using System;
using System.Collections.Generic;
using System.Linq;

namespace ptmrs
{
    public class PeptideSpectrumMatch
    {
        public static bool AllowNlWithEcdEtd;
        private int _charge;
        private double _precMassZ;
        internal List<FragmentIonType> _baseTypes;
        private List<Peak> _peaks;
        private List<FragmentIonType> _allFITs = new List<FragmentIonType>();
        private AminoAcidSequence _sequence;

        public string SequenceOneLetterCode
        {
            get
            {
                if (_sequence != null)
                    return _sequence.OneLetterCodeString();
                return "";
            }
        }

        internal FitComposition FragmentIonTypeComposition
        {
            get
            {
                switch (SpectrumType)
                {
                    case SpectrumType.None:
                        throw new ArgumentException("Spectrum type no set");
                    case SpectrumType.CID_CAD:
                    case SpectrumType.HCD:
                        return FitComposition.B | FitComposition.Y;
                    case SpectrumType.ECD_ETD:
                        return FitComposition.C | FitComposition.ZRadical | FitComposition.ZPrime;
                    case SpectrumType.EThcD:
                        return FitComposition.B | FitComposition.Y | FitComposition.C | FitComposition.ZRadical |
                               FitComposition.ZPrime;
                    default:
                        throw new ArgumentException("Spectrum type unknown");
                }
            }
        }

        internal int MaxFragmentIonCharge
        {
            get
            {
                switch (SpectrumType)
                {
                    case SpectrumType.None:
                        throw new ArgumentException("Spectrum type no set");
                    case SpectrumType.CID_CAD:
                    case SpectrumType.HCD:
                    case SpectrumType.EThcD:
                        return Math.Min(Charge, 2);
                    case SpectrumType.ECD_ETD:
                        return 1;
                    default:
                        throw new ArgumentException("Spectrum type unknown");
                }
            }
        }

        public int Charge => _charge;

        public SpectrumType SpectrumType { get; set; }

        /// <summary>
        /// Precursor m/z
        /// </summary>
        public double PrecMassZ => _precMassZ;

        public int Id { get; set; }

        public List<Peak> Peaks
        {
            get { return _peaks; }
            set
            {
                _peaks = value;
                _peaks.Sort();
            }
        }

        public AminoAcidSequence Sequence
        {
            get { return _sequence; }
            private set
            {
                if (value == null) throw new ArgumentNullException(nameof(value));
                _sequence = value;
            }
        }

        public PeptideSpectrumMatch(int id, SpectrumType type, int charge, double precMassZ, Peak[] peaks,
            AminoAcidSequence sequence)
        {
            Id = id;
            _charge = charge;
            _precMassZ = precMassZ;
            _sequence = sequence;
            _peaks = peaks.ToList();
            _peaks.Sort();
            SpectrumType = type;
        }

        public static SpectrumType parseType(string str)
        {
            if (str.ToUpper().Equals("CID") || str.ToUpper().Equals("CAD"))
                return SpectrumType.CID_CAD;
            if (str.ToUpper().Equals("ECD") || str.ToUpper().Equals("ETD"))
                return SpectrumType.ECD_ETD;
            if (str.ToUpper().Equals("HCD"))
                return SpectrumType.HCD;
            if (str.ToUpper().Equals("ETHCD"))
                return SpectrumType.EThcD;
            throw new ArgumentException("Fragmentation method: \"" + str +
                                        "\" unequal to \"CID\", \"CAD\", \"ECD\", \"ETD\" and \"HCD\"");
        }

        private static bool AddFITtoList(FragmentIonType fit, List<FragmentIonType> list)
        {
            if (fit == null)
            {
                return false;
            }
            if (list == null)
            {
                list = new List<FragmentIonType>();
            }
            for (int i = 0; i < list.Count; i++)
            {
                if (list[i].Equals(fit, true, false, false, false, false))
                {
                    return false;
                }
            }
            list.Add(fit);
            return true;
        }

        private static bool AddFITtoList(List<FragmentIonType> fits, List<FragmentIonType> list)
        {
            if (fits == null || fits.Count <= 0)
            {
                return false;
            }
            if (list == null)
            {
                list = new List<FragmentIonType>();
            }
            foreach (FragmentIonType current in fits)
            {
                for (int i = 0; i < list.Count; i++)
                {
                    list[i].Equals(current, true, false, false, false, false);
                }
                list.Add(current);
            }
            return true;
        }

        internal List<FragmentIonType> GetBasisFITs(FitComposition usedFITs)
        {
            int maxFragmentIonCharge = MaxFragmentIonCharge;
            List<FragmentIonType> list = new List<FragmentIonType>(maxFragmentIonCharge * usedFITs.Count());
            for (int i = 1; i <= maxFragmentIonCharge; i++)
            {
                if (usedFITs.HasFlag(FitComposition.B))
                    list.Add(GetBFragmentIonType(i));
                if (usedFITs.HasFlag(FitComposition.Y))
                    list.Add(GetYFragmentIonType(i));
            }
            if (usedFITs.HasFlag(FitComposition.C))
                list.Add(GetCFragmentIonType(1));
            if (usedFITs.HasFlag(FitComposition.ZPrime))
                list.Add(GetZPrimeFragmentIonType(1));
            if (usedFITs.HasFlag(FitComposition.ZRadical))
                list.Add(GetZRadicalFragmentIonType(1));
            return list;
        }

        internal static List<FragmentIonType> GetAllFiTs(AminoAcidSequence sequence,
            List<AminoAcidModification> scoredAaMods, List<FragmentIonType> baseFiTs, FitComposition considerNlPeaks,
            bool consideredOnly)
        {
            if (baseFiTs == null)
                return null;
            if (considerNlPeaks == FitComposition.No)
                return baseFiTs;

            List<FragmentIonType> list = new List<FragmentIonType>();
            Dictionary<int, Tuple<int, AminoAcidModification>> modificationAaSeqPosPairs =
                sequence.GetModificationAASeqPosPairs(true, scoredAaMods);
            if (modificationAaSeqPosPairs == null)
                return baseFiTs;

            List<int> positions = modificationAaSeqPosPairs.Select(w => w.Key).ToList();
            List<AminoAcidModification> mods = modificationAaSeqPosPairs.Select(w => w.Value.Item2).ToList();
            List<int> list3 = modificationAaSeqPosPairs.Select(w => w.Value.Item1).ToList();

            List<Tuple<int, AminoAcidModification.NeutralLoss>> list4 =
                new List<Tuple<int, AminoAcidModification.NeutralLoss>>();
            int i;
            for (i = 0; i < mods.Count; i++)
            {
                AminoAcid aminoAcid = sequence.GetAminoAcid( <>
                c__DisplayClassb.positions[i])
                ;
                list4.AddRange(from w in mods[i].NeutralLosses(list3[i], aminoAcid.OneLetterCode)
                    select new Tuple<int, AminoAcidModification.NeutralLoss>( < > c__DisplayClassb.positions[i],
                w))
                ;
                foreach (FragmentIonType current in from w in baseFiTs
                    where  <>

                c__DisplayClassb.considerNLPeaks.HasFlag(w.FITComposition)
                select w)
                {
                    PeptideSpectrumMatch.AddFITtoList(
                            mods[i].getNLFragmentIonTypes(aminoAcid.OneLetterCode, list3[i], current,  < >
                        c__DisplayClassb.positions[i]),
                    list)
                    ;
                    if (!consideredOnly)
                    {
                        PeptideSpectrumMatch.AddFITtoList(
                                mods[i].getNLFragmentIonTypes(aminoAcid.OneLetterCode, list3[i], current,  < >
                            c__DisplayClassb.positions[i],
                        -18.0105646836,
                        "-H2O"),
                        list)
                        ;
                    }
                }
            }
            List<string> list5 = new List<string>();
            int num = 0;
            int num2 = 2147483647;
            for (int j = 2; j <= list4.Count; j++)
            {
                List<List<Tuple<int, AminoAcidModification.NeutralLoss>>> combinationsWithoutReplication =
                    list4.GetCombinationsWithoutReplication(j);
                foreach (List<Tuple<int, AminoAcidModification.NeutralLoss>> current2 in combinationsWithoutReplication)
                {
                    double num3 = 0.0;
                    num = 0;
                    num2 = 2147483647;
                    list5.Clear();
                    foreach (Tuple<int, AminoAcidModification.NeutralLoss> current3 in current2)
                    {
                        num = Math.Max(current3.Item1, num);
                        num2 = Math.Min(current3.Item1, num2);
                        num3 += current3.Item2.DeltaMass;
                        list5.Add(current3.Item2.Abbreviation);
                    }
                    string text = AminoAcidModification.FormatAbbrList(list5);
                    foreach (FragmentIonType current4 in from w in baseFiTs
                        where  <>

                    c__DisplayClassb.considerNLPeaks.HasFlag(w.FITComposition)
                    select w)
                    {
                        PeptideSpectrumMatch.AddFITtoList(
                            new FragmentIonType(current4, true, current4.IsNTerminal ? num : num2, -num3, text, "-"),
                            list);
                        if (!consideredOnly)
                        {
                            PeptideSpectrumMatch.AddFITtoList(
                                new FragmentIonType(current4, true, current4.IsNTerminal ? num : num2,
                                    -num3 - 18.0105646836, text + "-H2O", "-"), list);
                        }
                    }
                }
            }
            List<FragmentIonType> list6 = new List<FragmentIonType>(list.Count + baseFiTs.Count);
            foreach (FragmentIonType current5 in list)
            {
                int num4 = list6.TryGetFITIndex(current5);
                if (num4 < 0)
                {
                    list6.Add(current5);
                }
                else
                {
                    FragmentIonType fragmentIonType = list6[num4];
                    if ((current5.IsNTerminal && current5.Expansion < fragmentIonType.Expansion) ||
                        (current5.IsCTerminal && current5.Expansion > fragmentIonType.Expansion))
                    {
                        list6[num4] = current5;
                    }
                }
            }
            list6.AddRange(baseFiTs);
            list6.TrimExcess();
            return list6;
        }

        internal FragmentIonType GetBFragmentIonType(int z)
        {
            return new FragmentIonType(FragmentIonTerminalType.NTerminal, false, 0, z, 0-Atoms.HydrogenMass, FitComposition.B,
                string.Empty);
        }

        internal FragmentIonType GetYFragmentIonType(int z)
        {
            return new FragmentIonType(FragmentIonTerminalType.CTerminal, false, this._sequence.Length - 1, z,
                1.007825032, FitComposition.Y, string.Empty);
        }

        internal FragmentIonType GetCFragmentIonType(int z)
        {
            return new FragmentIonType(FragmentIonTerminalType.NTerminal, false, 0, z, 16.018724069, FitComposition.C,
                string.Empty);
        }

        internal FragmentIonType GetZRadicalFragmentIonType(int z)
        {
            return new FragmentIonType(FragmentIonTerminalType.CTerminal, false, _sequence.Length - 1, z,
                -15.010899037000001, FitComposition.ZRadical, string.Empty);
        }

        internal FragmentIonType GetZPrimeFragmentIonType(int z)
        {
            return new FragmentIonType(FragmentIonTerminalType.CTerminal, false, _sequence.Length - 1, z,
                -14.003074005000002, FitComposition.ZPrime, "+H");
        }

        public override bool Equals(object obj)
        {
            PeptideSpectrumMatch peptideSpectrumMatch = (PeptideSpectrumMatch) obj;
            return peptideSpectrumMatch != null &&
                   (SpectrumType == peptideSpectrumMatch.SpectrumType && _charge == peptideSpectrumMatch._charge &&
                    PtmMathHelper.Equal(_precMassZ, peptideSpectrumMatch.PrecMassZ, 1E-06) &&
                    _peaks.Equals(peptideSpectrumMatch.Peaks)) &&
                   _sequence.Equals(peptideSpectrumMatch.Sequence);
        }

        public override int GetHashCode()
        {
            double num = (double) SpectrumType.GetHashCode();
            num += _charge.GetHashCode() * 2;
            num += Math.Round(_precMassZ, 6).GetHashCode() * 2;
            num += _peaks.GetHashCode() * 2;
            num += _allFITs.GetHashCode() * 3;
            num += _sequence.GetHashCode() * 2;
            return (int) num % 2147483647;
        }

        // 4.20
        internal void ClearPeakAnnotations()
        {
            if (_peaks != null)
            {
                foreach (Peak current in _peaks)
                {
                    current.ClearMatchData();
                }
            }
        }

        internal static void MatchIsoformToPeaks(AminoAcidSequence potentialSequence, List<Peak> peaks,
            double massTolerance, SpectrumType specType, List<FragmentIonType> FITs)
        {
            int num = 0;
            List<FragmentIon> list = potentialSequence.CalcTheoreticalFragIons(specType, FITs);
            foreach (FragmentIon current in list)
            {
                num++;
                Peak.MatchPeak(potentialSequence, current, massTolerance, peaks);
            }
            potentialSequence.ConsiredTheoreticalPeakCount = num;
        }

        /// <summary>
        /// 查找肽段和谱峰的所有匹配离子。4-20
        /// </summary>
        /// <param name="potentialSequence">待查找肽段</param>
        /// <param name="peaks">待匹配谱峰</param>
        /// <param name="massTolerance">质量精度</param>
        /// <param name="specType">谱图类型</param>
        /// <param name="fiTs">碎片离子类型</param>
        /// <returns>所有匹配的离子及对应离子强度</returns>
        internal static List<Tuple<FragmentIon, double>> GetMatchingFragmentIons(AminoAcidSequence potentialSequence,
            List<Peak> peaks, double massTolerance, SpectrumType specType, List<FragmentIonType> fiTs)
        {
            List<Tuple<FragmentIon, double>> tupleList = new List<Tuple<FragmentIon, double>>();
            List<FragmentIon> theoreticalFragIons = potentialSequence.CalcTheoreticalFragIons(specType, fiTs);
            foreach (FragmentIon current in theoreticalFragIons)
            {
                double intensity = Peak.MatchPeak(current, massTolerance, peaks);
                if (intensity > 0)
                    tupleList.Add(new Tuple<FragmentIon, double>(current, intensity));
            }
            return tupleList;
        }

        internal static void InitializeMatchPeaksCount(AminoAcidSequence sequence, PeakExtractor peakExtractor,
            int minPeakDepth, int maxPeakDepth)
        {
            List<List<AminoAcidSequence.MatchData>> list = new List<List<AminoAcidSequence.MatchData>>(peakExtractor.MzWindowCount + 1);
            for (int i = 0; i < peakExtractor.MzWindowCount; i++)
            {
                int num = 0;
                int num2 = 0;
                List<AminoAcidSequence.MatchData> list2 =
                    new List<AminoAcidSequence.MatchData>(maxPeakDepth - minPeakDepth);
                for (int j = 1; j <= maxPeakDepth; j++)
                {
                    if (num + 1 < j)
                        list2.Add(new AminoAcidSequence.MatchData(num2, num));
                    else
                    {
                        Peak nMostIntensePeaks = peakExtractor.GetNMostIntensePeaks(i, j);
                        if (nMostIntensePeaks != null)
                        {
                            num++;
                            num2 += nMostIntensePeaks.IsMatchingAASequence(sequence);
                        }
                        if (j >= minPeakDepth)
                        {
                            list2.Add(new AminoAcidSequence.MatchData(num2, num));
                        }
                    }
                }
                list.Add(list2);
            }
            sequence.initializeMatchedPeaksCount(list);
        }
    }
}