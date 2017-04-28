using System;
using System.Collections.Generic;
using System.Globalization;
using System.Linq;
using System.Text;
using System.Xml.Linq;
using ptmrs.Data;

namespace ptmrs
{
    public class AminoAcidModification
    {
        private static HashSet<int> IDs = new HashSet<int>();

        private int hashcode;

        private Dictionary<TargetInformation, List<NeutralLoss>> _nl;

        public AminoAcidModification(int id, string name, string abbreviation, double massDelta,
            List<AminoAcid> targets,
            Dictionary<TargetInformation, List<NeutralLoss>> neutralLosses)
        {
            lock (IDs)
            {
                if (!IDs.Contains(id))
                {
                    IDs.Add(id);
                }
            }
            ID = id;
            Name = name;
            Abbreviation = abbreviation;
            MassDelta = massDelta;
            TargetAminoAcids = targets;
            _nl = neutralLosses;
            DiagnosticIons = null;
            CalculateHash();
        }

        public AminoAcidModification(int id, string name, string abbreviation, double massDelta,
            List<AminoAcid> targetAAs,
            Dictionary<TargetInformation, List<NeutralLoss>> neutralLosses,
            List<DiagnosticIon> diagnosticIons) : this(id, name, abbreviation, massDelta, targetAAs, neutralLosses)
        {
            DiagnosticIons = diagnosticIons;
        }

        public AminoAcidModification(int id, string name, string abbreviation, double massDelta,
            AminoAcidModification basedModification)
            : this(id, name, abbreviation, massDelta, basedModification.TargetAminoAcids, basedModification._nl)
        {
        }

        public bool HasTargetAAs => TargetAminoAcids != null && TargetAminoAcids.Count > 0;

        public int ID { get; }

        public string Name { get; }

        public string Abbreviation { get; }

        public double MassDelta { get; }

        public List<AminoAcid> TargetAminoAcids { get; }

        public List<DiagnosticIon> DiagnosticIons { get; }

        public long UnimodAccession { get; set; }

        private void CalculateHash()
        {
            double num;
            if (!string.IsNullOrEmpty(Name))
                num = Name.GetHashCode();
            else if (!string.IsNullOrEmpty(Abbreviation))
                num = Abbreviation.GetHashCode();
            else
                num = ID.GetHashCode();
            try
            {
                num += Math.Pow(MassDelta, 2.0);
                if (_nl != null && _nl.Count > 0)
                {
                    num += _nl.Sum(w => w.Key.GetHashCode() + w.Value.Sum(v => v.GetHashCode()));
                }
                num %= int.MaxValue;
                hashcode = (int) num;
            }
            catch (Exception)
            {
                hashcode = (int) (num % int.MaxValue);
            }
        }

        public List<NeutralLoss> NeutralLosses(TargetInformation aa)
        {
            if (_nl == null)
                return new List<NeutralLoss>();

            List<NeutralLoss> list = null;
            List<NeutralLoss> list2;
            if (!aa.WithoutMultiplicity)
                _nl.TryGetValue(aa, out list);

            _nl.TryGetValue(new TargetInformation(-1, aa.Residue), out list2);
            if (list2 == null && list == null)
            {
                return new List<NeutralLoss>();
            }
            if (list2 == null)
            {
                return list;
            }
            if (list == null)
            {
                return list2;
            }
            return list.Concat(list2).ToList();
        }

        public List<NeutralLoss> NeutralLosses()
        {
            return _nl?.SelectMany(w => w.Value).ToList();
        }

        public List<NeutralLoss> NeutralLosses(int multiplicity, char AA)
        {
            return NeutralLosses(new TargetInformation(multiplicity, AA));
        }

        public string ToStringNeutralLosses()
        {
            if (NeutralLossCount() <= 0)
                return "";
            var stringBuilder = new StringBuilder();
            foreach (var current in _nl)
                stringBuilder.Append(string.Join(",", from v in current.Value select v.ToString()));
            return stringBuilder.ToString();
        }

        public string ToStringNeutralLosses(char aminoAcid, int numberOfModifications)
        {
            if (NeutralLossCount() <= 0)
                return "";
            StringBuilder stringBuilder = new StringBuilder();
            foreach (var current in from w in _nl
                where w.Key.Residue == aminoAcid
                      && (w.Key.WithoutMultiplicity || w.Key.NumberOfModifications == numberOfModifications)
                select w)
                stringBuilder.Append(string.Join(",", current.Value.Select(v => v.ToString())));
            return stringBuilder.ToString();
        }

        public int NeutralLossCount()
        {
            if (_nl != null)
                return _nl.Sum(w => w.Value.Count);
            return 0;
        }

        public int NeutralLossCount(char aminoAcid, int numberOfModifications)
        {
            if (_nl == null)
                return 0;
            return (from w in _nl
                where w.Key.Residue == aminoAcid &&
                      (w.Key.WithoutMultiplicity || w.Key.NumberOfModifications == numberOfModifications)
                select w).Sum(w => w.Value.Count);
        }

        public int NeutralLossCount(TargetInformation info)
        {
            if (_nl == null)
                return 0;
            var num = 0;
            List<NeutralLoss> list;
            if (_nl.TryGetValue(info, out list))
                num = list.Count;
            List<NeutralLoss> list2;
            if (_nl.TryGetValue(new TargetInformation(-1, info.Residue), out list2))
                num += list2.Count;
            return num;
        }

        public void AddNeutralLosses(TargetInformation information, NeutralLoss nl)
        {
            if (_nl == null)
            {
                _nl = new Dictionary<TargetInformation, List<NeutralLoss>>();
            }
            List<NeutralLoss> list;
            if (_nl.TryGetValue(information, out list))
            {
                if (list.All(w => !PtmMathHelper.Equal(w.DeltaMass, nl.DeltaMass, 0.0001)))
                    list.Add(nl);
            }
            else if (!information.WithoutMultiplicity &&
                     _nl.TryGetValue(new TargetInformation(-1, information.Residue), out list))
            {
                if (list.All(w => !PtmMathHelper.Equal(w.DeltaMass, nl.DeltaMass, 0.0001)))
                    list.Add(nl);
            }
            else
                _nl.Add(information, new List<NeutralLoss>(1) {nl});
            CalculateHash();
        }

        public void ClearNeutralLosses()
        {
            _nl = null;
            CalculateHash();
        }

        public IEnumerator<KeyValuePair<TargetInformation, List<NeutralLoss>>> GetNeutalLossEnumerator()
        {
            if (_nl != null)
            {
                foreach (var current in _nl)
                {
                    yield return current;
                }
            }
        }

        public static int getUniqueID()
        {
            lock (IDs)
            {
                var num = int.MaxValue;
                while (IDs.Contains(num))
                {
                    num--;
                }
                return num;
            }
        }

        public override int GetHashCode()
        {
            return hashcode;
        }

        public override bool Equals(object obj)
        {
            try
            {
                if (obj == null || obj.GetType() != typeof(AminoAcidModification))
                    return false;

                return Equals((AminoAcidModification) obj);
            }
            catch (Exception)
            {
                return false;
            }
        }

        public bool Equals(AminoAcidModification value)
        {
            try
            {
                if (value == null || value.GetHashCode() != GetHashCode())
                {
                    return false;
                }
                List<AminoAcid> targets = value.TargetAminoAcids;
                return TargetAminoAcids.Any(acid => targets.Contains(acid));
            }
            catch (Exception)
            {
                return false;
            }
        }

        public override string ToString()
        {
            return Name;
        }

        public bool AllowsNL(char aa, int multiplicity)
        {
            return NeutralLossCount(new TargetInformation(multiplicity, aa)) > 0;
        }

        public AminoAcidModification Clone()
        {
            return new AminoAcidModification(ID, Name, Abbreviation, MassDelta, TargetAminoAcids, _nl);
        }

        public static bool FindAAMod(char id, List<AminoAcidModification> allMods, out AminoAcidModification outResult)
        {
            for (var i = 0; i < allMods.Count; i++)
            {
                var aminoAcidModification = allMods[i];
                if (id == aminoAcidModification.ID)
                {
                    outResult = aminoAcidModification.Clone();
                    return true;
                }
            }
            outResult = null;
            return false;
        }

        public static string FormatAbbrList(List<string> abbreviations)
        {
            StringBuilder stringBuilder = new StringBuilder();
            for (int index = 0; index < abbreviations.Count; index++)
            {
                String str = abbreviations[index];
                if (!string.IsNullOrEmpty(str))
                {
                    var num = (from w in abbreviations
                        where w.Equals(str)
                        select w).Count();
                    stringBuilder.Append(" -");
                    if (num > 1)
                    {
                        stringBuilder.Append(num);
                    }
                    stringBuilder.Append(str);
                }
            }
            return stringBuilder.ToString();
        }

        public void AddTargetAminoAcids(List<AminoAcid> targetAAs)
        {
            TargetAminoAcids.AddRange(targetAAs);
        }

        internal List<FragmentIonType> getNLFragmentIonTypes(char aa, int multiplicity, FragmentIonType baseFIT,
            int extension)
        {
            return getNLFragmentIonTypes(new TargetInformation(multiplicity, aa), baseFIT, extension);
        }

        internal List<FragmentIonType> getNLFragmentIonTypes(TargetInformation AA,
            FragmentIonType baseFIT, int Extension)
        {
            List<NeutralLoss> source = NeutralLosses(AA);
            if (source.Count == 0)
                return new List<FragmentIonType>(0);
            return source.Select(w => w.getFragmentIonType(baseFIT, Extension)).ToList();
        }

        internal List<FragmentIonType> getNLFragmentIonTypes(TargetInformation AA,
            FragmentIonType baseFIT, int Extension, double additionalDelta, string additionalAppentix)
        {
            List<NeutralLoss> source = NeutralLosses(AA);
            if (source.Count == 0)
                return new List<FragmentIonType>(0);

            return source.Select(w => w.getFragmentIonType(baseFIT, Extension, additionalDelta, additionalAppentix))
                .ToList();
        }

        internal List<FragmentIonType> getNLFragmentIonTypes(char aa, int multiplicity, FragmentIonType baseFIT,
            int extension, double additionalDelta, string additionalAppentix)
        {
            var nLFragmentIonTypes = getNLFragmentIonTypes(new TargetInformation(multiplicity, aa), baseFIT,
                extension, additionalDelta, additionalAppentix);
            return nLFragmentIonTypes;
        }

        // 2017-04-19
        public class NeutralLoss
        {
            public NeutralLoss(string abbreviation, double deltaMass)
            {
                DeltaMass = deltaMass;
                Abbreviation = abbreviation;
            }

            public NeutralLoss(double deltaMass) : this(null, deltaMass)
            {
            }

            public NeutralLoss(XElement neutralLoss) : this(
                neutralLoss.Attribute("abbreviation").Value, Convert.ToDouble(neutralLoss.Attribute("mass").Value))
            {
            }

            public double DeltaMass { get; }

            public string Abbreviation { get; }

            internal FragmentIonType getFragmentIonType(FragmentIonType baseFIT, int Extension)
            {
                return new FragmentIonType(baseFIT, Extension, -DeltaMass, Abbreviation);
            }

            internal FragmentIonType getFragmentIonType(FragmentIonType baseFIT, int Extension,
                double additionalMassDiff, string additionalAppentix)
            {
                if (string.IsNullOrEmpty(additionalAppentix))
                {
                    return new FragmentIonType(baseFIT, Extension, -DeltaMass + additionalMassDiff, Abbreviation);
                }
                return new FragmentIonType(baseFIT, Extension, -DeltaMass + additionalMassDiff,
                    Abbreviation + additionalAppentix);
            }

            public override string ToString()
            {
                if (string.IsNullOrEmpty(Abbreviation))
                {
                    return DeltaMass.ToString(CultureInfo.InvariantCulture);
                }
                return string.Join(" ", Abbreviation, "-", DeltaMass.ToString(CultureInfo.InvariantCulture));
            }

            public override bool Equals(object obj)
            {
                if (obj == null || GetType() != obj.GetType())
                {
                    return false;
                }
                var neutralLoss = obj as NeutralLoss;
                return Math.Abs(DeltaMass - neutralLoss.DeltaMass) <= 1E-06 &&
                       string.CompareOrdinal(neutralLoss.Abbreviation, Abbreviation) == 0;
            }

            public override int GetHashCode()
            {
                return (int) DeltaMass * 1000;
            }
        }

        // 2017-04-19
        public class TargetInformation
        {
            public const int ForAllNumberOfModifications = -1;
            private int _numberOfModifications; // 对二甲基，其值为2

            public TargetInformation(int multiplicity, char residue)
            {
                _numberOfModifications = multiplicity;
                Residue = residue;
            }

            public TargetInformation(char residue) : this(-1, residue)
            {
            }

            /// <summary>
            /// 修改的氨基酸位点
            /// </summary>
            public char Residue { get; }

            public bool WithoutMultiplicity => _numberOfModifications <= 0;

            /// <summary>
            /// 多二甲基，其值为2
            /// </summary>
            public int NumberOfModifications
            {
                get
                {
                    if (_numberOfModifications <= 0)
                    {
                        return 0;
                    }
                    return _numberOfModifications;
                }
            }

            /// <summary>
            /// 解析配置文件
            /// </summary>
            /// <param name="targetInformationXML"></param>
            /// <param name="Targets"></param>
            /// <returns></returns>
            public static List<TargetInformation> CreateNew(XElement targetInformationXML, List<AminoAcid> Targets)
            {
                // 没有任何氨基酸位点
                if (targetInformationXML.Elements("Target").Count() == 0)
                {
                    return (from aa in Targets
                        select new TargetInformation(-1, aa.OneLetterCode)).ToList();
                }
                var list = new List<TargetInformation>();
                foreach (XElement current in targetInformationXML.Elements("Target"))
                {
                    char c = Convert.ToChar(current.Attribute("aminoacid").Value, Constants.EN_US_CULTURE_INFO);
                    var multiplicity = -1;
                    if (current.Attribute("factor") != null &&
                        !int.TryParse(current.Attribute("factor").Value, out multiplicity))
                    {
                        multiplicity = -1;
                    }
                    if (AminoAcid.GetAminoAcid(c) != null)
                    {
                        list.Add(new TargetInformation(multiplicity, c));
                    }
                }
                return list;
            }

            public override bool Equals(object obj)
            {
                if (obj == null || GetType() != obj.GetType())
                {
                    return false;
                }
                var targetInformation = obj as TargetInformation;
                return _numberOfModifications == targetInformation._numberOfModifications &&
                       Residue == targetInformation.Residue;
            }

            public override int GetHashCode()
            {
                return (int) (_numberOfModifications * (double) Residue.GetHashCode());
            }

            public override string ToString()
            {
                if (_numberOfModifications <= 0)
                {
                    return $"on {Residue}";
                }
                return $"x{_numberOfModifications} on {Residue}";
            }
        }
    }
}