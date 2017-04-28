using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace ptmrs
{
    public class AminoAcidSequence
    {
        // 4.19
        public struct MatchData
        {
            public int MatchedPeakCount;
            public int ExtractedPeakCount;

            /// <summary>
            /// Constructor.
            /// </summary>
            /// <param name="matchedPeakCount">匹配的谱峰数目</param>
            /// <param name="extractedPeakCount">提取的谱峰数目</param>
            public MatchData(int matchedPeakCount, int extractedPeakCount)
            {
                MatchedPeakCount = matchedPeakCount;
                ExtractedPeakCount = extractedPeakCount;
            }
        }

        private readonly int _id;
        private readonly Terminus _nTerminus;
        private readonly Terminus _cTerminus;

        private readonly List<AminoAcid> _AminoAcids;

        // 保存该肽段在所有窗口的匹配谱峰信息 windowid -> peakdepth-id
        private List<List<MatchData>> _matchedPeaksCount;

        private string _myString;

        public int ConsiredTheoreticalPeakCount { get; set; }
        public int Rank { get; set; }
        public int ID => _id;

        public Terminus NTerminus => _nTerminus;
        public Terminus CTerminus => _cTerminus;

        /// <summary>
        /// 肽段长度，或氨基酸数目
        /// </summary>
        public int Length => _AminoAcids.Count;

        /// <summary>
        /// 获取该肽段在指定窗口指定 peak depth 下的匹配谱峰和提取谱峰数目  2017-4-20
        /// </summary>
        /// <param name="window">窗口 index</param>
        /// <param name="peakDepth">peak depth</param>
        /// <param name="minPeakDepth">最小 peak depth 值</param>
        /// <returns>肽段在指定窗口指定 peakdepth 下的匹配谱峰信息</returns>
        public MatchData getMatchData(int window, int peakDepth, int minPeakDepth)
        {
            return _matchedPeaksCount[window][peakDepth - minPeakDepth];
        }

        /// <summary>
        /// 统计该肽段在所有窗口匹配和提取的谱峰数目 // 2017-4-20
        /// </summary>
        /// <param name="peakDepth">保存该肽段在不同窗口的 peakdepth</param>
        /// <param name="minPeakDepth"></param>
        /// <returns></returns>
        public MatchData getMatchData(List<int> peakDepth, int minPeakDepth)
        {
            MatchData matchData = new MatchData();
            for (int i = 0; i < peakDepth.Count; i++)
            {
                MatchData matchData2 = _matchedPeaksCount[i][peakDepth[i] - minPeakDepth];
                matchData.ExtractedPeakCount += matchData2.ExtractedPeakCount;
                matchData.MatchedPeakCount += matchData2.MatchedPeakCount;
            }
            return matchData;
        }

        /// <summary>
        /// 清除所有匹配谱峰信息.
        /// </summary>
        public void ClearMatchData()
        {
            ConsiredTheoreticalPeakCount = 0;
            _matchedPeaksCount = null;
        }

        internal static void MatchIsoformToPeaks(AminoAcidSequence potentialSequence, List<Peak> peaks,
            double massTolerance, SpectrumType specType, List<FragmentIonType> fiTs)
        {
            int num = 0;
            List<FragmentIon> list = potentialSequence.CalcTheoreticalFragIons(specType, fiTs);
            foreach (FragmentIon current in list)
            {
                num++;
                Peak.MatchPeak(potentialSequence, current, massTolerance, peaks);
            }
            potentialSequence.ConsiredTheoreticalPeakCount = num;
        }

        /// <summary>
        /// 初始化匹配的数据
        /// </summary>
        public void initializeMatchedPeaksCount(List<List<MatchData>> matchedPeaksCount)
        {
            if (_matchedPeaksCount != null)
            {
                throw new ArgumentException("Match peak count dict already initilized");
            }
            _matchedPeaksCount = matchedPeaksCount;
        }

        // 2017-4-20
        public AminoAcidSequence() : this(-1, -1, new List<AminoAcid>(), new Terminus(), new Terminus())
        {
        }

        // 2017-4-20
        public AminoAcidSequence(int id, int rank, List<AminoAcid> aaSequence, Terminus n, Terminus c)
        {
            _id = id;
            _nTerminus = n;
            _cTerminus = c;
            _AminoAcids = aaSequence;
            Rank = rank;
        }

        public static AminoAcidSequence Create(int ID, int rank, string sequence, List<AminoAcidModification> mods,
            string modPos)
        {
            List<AminoAcid> list = ParseAASequence(sequence);
            Terminus terminus = new Terminus(FragmentIonTerminalType.NTerminal);
            char c = modPos[0];
            if (c != '0')
            {
                AminoAcidModification modification;
                if (AminoAcidModification.FindAAMod(c, mods, out modification))
                    terminus.setModification(modification, 1);
            }
            Terminus terminus2 = new Terminus(FragmentIonTerminalType.CTerminal);
            c = modPos[modPos.Length - 1];
            if (c != '0')
            {
                AminoAcidModification modification;
                if (AminoAcidModification.FindAAMod(c, mods, out modification))
                    terminus2.setModification(modification, 1);
            }
            string text = modPos.Substring(modPos.IndexOf('.') + 1, modPos.LastIndexOf('.'));
            for (int i = 0; i < text.Length; i++)
            {
                c = text[i];
                if (c != '0')
                {
                    AminoAcidModification modification;
                    if (AminoAcidModification.FindAAMod(c, mods, out modification))
                    {
                        AminoAcid aminoAcid = list[i];
                        if (aminoAcid != null)
                            aminoAcid.setModification(modification, 1);
                    }
                }
            }
            return new AminoAcidSequence(ID, rank, list, terminus, terminus2);
        }

        public override string ToString()
        {
            if (this._myString == null)
            {
                StringBuilder stringBuilder = new StringBuilder();
                stringBuilder.Append("ID: ");
                stringBuilder.Append(this._id);
                stringBuilder.Append(" ");
                if (this.NTerminus != null)
                {
                    stringBuilder.Append(this.NTerminus);
                }
                else
                {
                    stringBuilder.Append("N");
                }
                stringBuilder.Append("-");
                foreach (AminoAcid current in this._AminoAcids)
                {
                    stringBuilder.Append(current);
                }
                stringBuilder.Append("-");
                if (this.CTerminus != null)
                {
                    stringBuilder.Append(this.CTerminus);
                }
                else
                {
                    stringBuilder.Append("C");
                }
                stringBuilder.Append(" ");
                foreach (var current2 in from aa in _AminoAcids.Select((AminoAcid amino, int index) => new
                    {
                        amino,
                        index
                    })
                    where aa.amino.IsModified
                    select aa)
                {
                    stringBuilder.Append(current2.index + 1);
                    stringBuilder.Append(current2.amino.GetModificationString());
                    stringBuilder.Append(" ");
                }
                this._myString = stringBuilder.ToString();
            }
            return this._myString;
        }

        public override bool Equals(object obj)
        {
            AminoAcidSequence aminoAcidSequence = (AminoAcidSequence) obj;
            return aminoAcidSequence != null && this._id == aminoAcidSequence._id;
        }

        public override int GetHashCode()
        {
            return this._id;
        }

        public bool IsEmpty()
        {
            return this._AminoAcids.Count == 0;
        }

        // 4.19
        public double CalcTheorPrecursorMr()
        {
            return _nTerminus.GetTotalMass() + AminoAcid.SummAAMasses(_AminoAcids) + _cTerminus.GetTotalMass();
        }

        public static List<AminoAcid> ParseAASequence(string sequence)
        {
            List<AminoAcid> list = new List<AminoAcid>(sequence.Length);
            for (int i = 0; i < sequence.Length; i++)
            {
                char letter = sequence[i];
                AminoAcid item;
                if (AminoAcid.FindAAResidue(letter, out item))
                {
                    list.Add(item);
                }
            }
            list.TrimExcess();
            return list;
        }

        public static List<AminoAcid> ParseAASequence(List<char> sequence)
        {
            List<AminoAcid> list = new List<AminoAcid>(sequence.Count);
            foreach (char current in sequence)
            {
                AminoAcid aminoAcid = AminoAcid.GetAminoAcid(current);
                if (aminoAcid == null)
                {
                    throw new ArgumentNullException("AminoAcid not found. Please contact your Administrator");
                }
                list.Add(aminoAcid);
            }
            list.TrimExcess();
            return list;
        }

        public string OneLetterCodeString()
        {
            return string.Join<char>("", (from w in this._AminoAcids
                select w.OneLetterCode).ToArray<char>());
        }

        public List<int> GetPositionsOf(List<AminoAcid> targetAminoAcids, bool ignoreCTerminalKandR, bool basedOnZero)
        {
            if (_AminoAcids == null)
            {
                return null;
            }
            List<int> list = new List<int>();
            if (basedOnZero)
            {
                for (int i = 0; i < this._AminoAcids.Count; i++)
                {
                    foreach (AminoAcid current in targetAminoAcids)
                    {
                        if (current.Equals(this._AminoAcids[i], true) &&
                            (!ignoreCTerminalKandR || i != this._AminoAcids.Count - 1 ||
                             (this._AminoAcids[i].OneLetterCode != 'K' && this._AminoAcids[i].OneLetterCode != 'R')))
                        {
                            list.Add(i);
                            break;
                        }
                    }
                }
            }
            else
            {
                for (int j = 0; j < _AminoAcids.Count; j++)
                {
                    foreach (AminoAcid current2 in targetAminoAcids)
                    {
                        if (current2.Equals(_AminoAcids[j], true) &&
                            (!ignoreCTerminalKandR || j != _AminoAcids.Count - 1 ||
                             (_AminoAcids[j].OneLetterCode != 'K' && _AminoAcids[j].OneLetterCode != 'R')))
                        {
                            list.Add(j + 1);
                            break;
                        }
                    }
                }
            }
            return list;
        }

        internal List<ModificationPosition> GetPositionsAndNumberOf(List<AminoAcid> targetAminoAcids,
            bool ignoreCTerminalKandR, bool basedOnZero)
        {
            if (_AminoAcids == null)
                return null;
            List<ModificationPosition> list = new List<ModificationPosition>();
            if (basedOnZero)
            {
                for (int i = 0; i < _AminoAcids.Count; i++)
                {
                    foreach (AminoAcid current in targetAminoAcids)
                    {
                        if (current.Equals(_AminoAcids[i], true) &&
                            (!ignoreCTerminalKandR || i != _AminoAcids.Count - 1 ||
                             (_AminoAcids[i].OneLetterCode != 'K' && _AminoAcids[i].OneLetterCode != 'R')))
                        {
                            list.Add(new ModificationPosition(i, this._AminoAcids[i].NumberOfModifications, true));
                            break;
                        }
                    }
                }
            }
            else
            {
                for (int j = 0; j < _AminoAcids.Count; j++)
                {
                    foreach (AminoAcid current2 in targetAminoAcids)
                    {
                        if (current2.Equals(_AminoAcids[j], true) &&
                            (!ignoreCTerminalKandR || j != _AminoAcids.Count - 1 ||
                             (_AminoAcids[j].OneLetterCode != 'K' && _AminoAcids[j].OneLetterCode != 'R')))
                        {
                            list.Add(new ModificationPosition(j + 1, _AminoAcids[j].NumberOfModifications,
                                false));
                            break;
                        }
                    }
                }
            }
            return list;
        }

        public int GetNumberSites(AminoAcidModification mod)
        {
            int result;
            int num = 0;
            using (IEnumerator<AminoAcid> enumerator = (from w in _AminoAcids
                where w.IsModified
                select w).GetEnumerator())
            {
                while (enumerator.MoveNext())
                {
                    AminoAcid aa = enumerator.Current;
                    if (aa.Modification.Equals(mod))
                    {
                        if (mod.TargetAminoAcids.Find(acid => acid.OneLetterCode == aa.OneLetterCode) !=
                            null)
                        {
                            num += aa.NumberOfModifications;
                        }
                    }
                }
            }
            result = num;
            return result;
        }

        internal List<ModificationPosition> GetPositionsOf(AminoAcidModification mod, bool basedOnZero)
        {
            List<ModificationPosition> list = new List<ModificationPosition>();
            if (basedOnZero)
            {
                int i;
                for (i = 0; i < _AminoAcids.Count; i++)
                {
                    if (_AminoAcids[i].IsModified && _AminoAcids[i].Modification.Equals(mod))
                    {
                        if (mod.TargetAminoAcids.Find(acid => acid.OneLetterCode == _AminoAcids[i].OneLetterCode) !=
                            null)
                            list.Add(new ModificationPosition(i, _AminoAcids[i].NumberOfModifications, true));
                    }
                }
            }
            else
            {
                int i;
                for (i = 0; i < _AminoAcids.Count; i++)
                {
                    if (_AminoAcids[i].IsModified && _AminoAcids[i].Modification.Equals(mod))
                    {
                        if (mod.TargetAminoAcids.Find(acid => acid.OneLetterCode == _AminoAcids[i].OneLetterCode) !=
                            null)
                            list.Add(new ModificationPosition(i + 1, _AminoAcids[i].NumberOfModifications, false));
                    }
                }
            }
            return list;
        }

        public AminoAcid GetAminoAcid(int position)
        {
            return _AminoAcids[position];
        }

        public List<AminoAcid> GetAminoAcids()
        {
            return _AminoAcids;
        }

        public Dictionary<int, Tuple<int, AminoAcidModification>> GetPositionsOf(List<AminoAcidModification> mods,
            bool basedOnZero)
        {
            if (_AminoAcids == null)
            {
                return new Dictionary<int, Tuple<int, AminoAcidModification>>();
            }
            return (from aa in _AminoAcids.Select((AminoAcid value, int index) => new
                {
                    value,
                    index
                })
                where aa.value.IsModified && mods.Any(
                          (AminoAcidModification aamod) => aamod.Equals(aa.value.Modification) &&
                                                           aamod.TargetAminoAcids.Contains(aa.value))
                select aa).ToDictionary(w => w.index,
                w => new Tuple<int, AminoAcidModification>(w.value.NumberOfModifications, w.value.Modification));
        }

        internal bool HasModification(ModificationPosition modificationPosition, AminoAcidModification mod)
        {
            return this._AminoAcids[modificationPosition.ZeroBasedPosition].IsModified &&
                   this._AminoAcids[modificationPosition.ZeroBasedPosition].NumberOfModifications ==
                   modificationPosition.NumberOfModifications &&
                   this._AminoAcids[modificationPosition.ZeroBasedPosition].Modification.Equals(mod) &&
                   mod.TargetAminoAcids.Contains(this._AminoAcids[modificationPosition.ZeroBasedPosition]);
        }

        public bool HasModification(int position, int atleastNModifications, AminoAcidModification mod)
        {
            return this._AminoAcids[position].IsModified &&
                   atleastNModifications <= this._AminoAcids[position].NumberOfModifications &&
                   this._AminoAcids[position].Modification.Equals(mod) &&
                   mod.TargetAminoAcids.Contains(this._AminoAcids[position]);
        }

        /// <summary>
        /// 如果包含特定修饰，返回 true.
        /// </summary>
        public bool HasModification(AminoAcidModification modification)
        {
            return _AminoAcids.Any(w => modification.Equals(w.Modification));
        }

        internal bool HasModification(PTMSiteProbability.SiteProbability siteProb)
        {
            bool result;
            if (siteProb.Position.OneBasedPosition > _AminoAcids.Count)
                result = false;
            else
            {
                AminoAcid aminoAcid = _AminoAcids[siteProb.Position.ZeroBasedPosition];
                result = aminoAcid.NumberOfModifications == siteProb.Position.NumberOfModifications &&
                         siteProb.Modification.Equals(aminoAcid.Modification);
            }
            return result;
        }

        public Dictionary<int, Tuple<int, AminoAcidModification>> GetModificationAASeqPosPairs(bool onlyWithNL)
        {
            return GetModificationAASeqPosPairs(onlyWithNL, null);
        }

        public Dictionary<int, Tuple<int, AminoAcidModification>> GetModificationAASeqPosPairs(bool onlyWithNL,
            List<AminoAcidModification> aaModType)
        {
            Dictionary<int, Tuple<int, AminoAcidModification>> dictionary =
                new Dictionary<int, Tuple<int, AminoAcidModification>>();

            foreach (AminoAcid aminoAcid in _AminoAcids)
            {
            }

            foreach (var current in from aaIndexPair in _AminoAcids.Select((value, index) => new
                {
                    aa = value,
                    pos = index
                })
                where aaIndexPair.aa.IsModified
                select aaIndexPair)
            {
                if (!onlyWithNL || current.aa.HasNL)
                {
                    AminoAcid aminoAcid = current.aa;
                    int pos = current.pos;
                    if (aaModType == null || aaModType.Any(aaMod =>
                    {
                        if (aminoAcid.Modification.Equals(aaMod))
                            return aaMod.TargetAminoAcids.Select(a => a.OneLetterCode)
                                .Contains(aminoAcid.OneLetterCode);
                        return false;
                    }))
                        dictionary.Add(pos,
                            new Tuple<int, AminoAcidModification>(aminoAcid.NumberOfModifications,
                                aminoAcid.Modification));
                }
            }
            return dictionary;
        }

        internal List<ModificationPosition> GetModificationMap(AminoAcidModification modification)
        {
            List<ModificationPosition> result;
            try
            {
                List<ModificationPosition> list = new List<ModificationPosition>(3);
                for (int i = 0; i < this._AminoAcids.Count; i++)
                {
                    if (modification.Equals(this._AminoAcids[i].Modification))
                    {
                        list.Add(new ModificationPosition(i, this._AminoAcids[i].NumberOfModifications, true));
                    }
                }
                list.TrimExcess();
                result = list;
            }
            catch (Exception)
            {
                throw;
            }
            return result;
        }

        public int Max_NumberOfModifications_perTarget(AminoAcidModification modification)
        {
            return (from w in this._AminoAcids
                where modification.Equals(w.Modification)
                select w).Max((AminoAcid w) => w.NumberOfModifications);
        }

        public double GetAverageTotalAminoAcidMass()
        {
            double num = 0.0;
            for (int i = 0; i < this._AminoAcids.Count; i++)
            {
                num += this._AminoAcids[i].GetTotalMass();
            }
            return num / (double) this._AminoAcids.Count;
        }

        public string GetUnmodifiedSequence()
        {
            StringBuilder stringBuilder = new StringBuilder();
            for (int i = 0; i < this._AminoAcids.Count; i++)
            {
                stringBuilder.Append(this._AminoAcids[i].OneLetterCode);
            }
            return stringBuilder.ToString();
        }

        public List<AminoAcid> GetUnmodifiedAASequence(List<AminoAcidModification> removeAAMod)
        {
            List<AminoAcid> list = new List<AminoAcid>(this._AminoAcids.Count);
            AminoAcidModification newModification;
            char oneletter;
            foreach (AminoAcid current in this._AminoAcids)
            {
                oneletter = current.OneLetterCode;
                newModification = current.Modification;
                int num = current.NumberOfModifications;
                if (num > 0)
                {
                    if (removeAAMod.Any(delegate(AminoAcidModification aaMod)
                    {
                        if (aaMod.Equals(newModification))
                        {
                            return (from w in aaMod.TargetAminoAcids
                                select w.OneLetterCode).Contains(oneletter);
                        }
                        return false;
                    }))
                    {
                        newModification = null;
                        num = 0;
                    }
                }
                list.Add(new AminoAcid(current.OneLetterCode, current.ThreeLetterCode, current.Name,
                    current.MonoisotopicMass, newModification, num));
            }
            return list;
        }

        public static List<AminoAcid> GetModifiedAASequence(List<AminoAcid> unmodifiedAASeq,
            Dictionary<int, Tuple<int, AminoAcidModification>> modifiedAASeqPos)
        {
            List<AminoAcid> list = new List<AminoAcid>(unmodifiedAASeq.Count);
            for (int i = 1; i <= unmodifiedAASeq.Count; i++)
            {
                AminoAcid aminoAcid = unmodifiedAASeq[i - 1];
                Tuple<int, AminoAcidModification> tuple;
                if (!modifiedAASeqPos.TryGetValue(i, out tuple))
                {
                    tuple =
                        new Tuple<int, AminoAcidModification>(aminoAcid.NumberOfModifications, aminoAcid.Modification);
                }
                list.Add(new AminoAcid(aminoAcid.OneLetterCode, aminoAcid.ThreeLetterCode, aminoAcid.Name,
                    aminoAcid.MonoisotopicMass, tuple.Item2, tuple.Item1));
            }
            return list;
        }

        /// <summary>
        /// 获得该肽段指定离子类型的所有碎片离子质量
        /// </summary>
        /// <param name="type">离子类型</param>
        /// <param name="avoidFragmentationNterminalPro">true 如果谱图类型为 ECD_ETD</param>
        /// <returns>该肽段的碎片离子质量</returns>
        internal FitMasses GetFitMasses(FragmentIonType type, bool avoidFragmentationNterminalPro)
        {
            int expansion = type.Expansion;
            if (type.IsNonTerminal)
            {
                double mass = (CalcTheorPrecursorMr() + type.MassDiff + Atoms.ProtonMass * type.Charge) / type.Charge;
                return new FitMasses(type, new List<FragmentIon>(1)
                {
                    new FragmentIon(mass, type, -1)
                });
            }
            FitMasses fItMasses = new FitMasses(type, new List<FragmentIon>(_AminoAcids.Count - 1));
            for (int i = 1; i < _AminoAcids.Count; i++)
            {
                if (type.IsNTerminal)
                {
                    if (expansion == FragmentIonType.FullExpansion || expansion + 1 <= i)
                    {
                        double mass;
                        if (!avoidFragmentationNterminalPro || _AminoAcids[i].OneLetterCode != 'P')
                        {
                            mass = (_nTerminus.GetTotalMass() + AminoAcid.SumAAMasses(i, _AminoAcids, false) +
                                    type.MassDiff + Atoms.ProtonMass * type.Charge) / type.Charge;
                        }
                        else
                            mass = -1.0;
                        fItMasses.FragmentIons.Add(new FragmentIon(mass, type, i));
                    }
                }
                else if (type.IsCTerminal && (expansion == -1 || expansion >= i))
                {
                    double mass;
                    if (!avoidFragmentationNterminalPro || !type.IsZIon() || _AminoAcids[i].OneLetterCode != 'P')
                    {
                        mass = (_cTerminus.GetTotalMass() +
                                AminoAcid.SumAAMasses(_AminoAcids.Count - i, _AminoAcids, true) +
                                type.MassDiff + Atoms.ProtonMass * type.Charge) / type.Charge;
                    }
                    else
                        mass = -1.0;
                    fItMasses.FragmentIons.Add(new FragmentIon(mass, type, _AminoAcids.Count - i));
                }
            }
            return fItMasses;
        }

        public static IEnumerable<List<AminoAcidSequence>> GetIsoformPeptideGroups(
            IEnumerable<AminoAcidSequence> peptides, List<Tuple<int, AminoAcidModification>> possibleModifications)
        {
            foreach (IGrouping<string, AminoAcidSequence> current in from g in peptides
                group g by g.OneLetterCodeString())
            {
                List<AminoAcidSequence> list = current.OrderBy(o => o.Rank).ToList();
                while (list.Count > 0)
                {
                    if (list.Count <= 1)
                    {
                        yield return list;
                        list.Clear();
                    }
                    else
                    {
                        Func<AminoAcidSequence, List<int>> getModificationIdsByAscendingOrder =
                            peptideItem => peptideItem.GetModificationAASeqPosPairs(false)
                                .SelectMany(s => Enumerable.Repeat(s.Value.Item2.ID, s.Value.Item1))
                                .OrderBy(o => o)
                                .ToList();
                        Func<AminoAcidSequence, int> getUnusedPTMModifications =
                            peptideItem => peptideItem.GetModificationAASeqPosPairs(false)
                                .Where(posAndMod => !
                                    possibleModifications.Any(delegate(Tuple<int, AminoAcidModification> w)
                                    {
                                        if (w.Item2.Equals(posAndMod.Value.Item2) &&
                                            w.Item1 == posAndMod.Value.Item1)
                                        {
                                            return w.Item2.TargetAminoAcids
                                                .Select(aa => aa.OneLetterCode)
                                                .Contains(peptideItem._AminoAcids[posAndMod.Key].OneLetterCode);
                                        }
                                        return false;
                                    }))
                                .Count();
                        List<int> modificationsOfCurrentIsoform =
                            getModificationIdsByAscendingOrder(list.First());
                        int numberOfUnusedMods = getUnusedPTMModifications(list.First());
                        List<AminoAcidSequence> list2 =
                            list.FindAll(r => getModificationIdsByAscendingOrder(r)
                                                  .SequenceEqual(modificationsOfCurrentIsoform) &&
                                              getUnusedPTMModifications(r) == numberOfUnusedMods);
                        list = list.Except(list2).ToList();
                        yield return list2;
                    }
                }
            }
            yield break;
        }

        internal List<FragmentIon> CalcTheoreticalFragIons(SpectrumType specType, List<FragmentIonType> fiTs)
        {
            return fiTs.SelectMany(fit => GetFitMasses(fit, specType == SpectrumType.ECD_ETD)
                    .FragmentIons.Where(fi => fi.Mass != -1.0))
                .ToList();
        }
    }
}