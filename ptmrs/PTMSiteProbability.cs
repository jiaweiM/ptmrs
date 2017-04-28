using System;
using System.Collections.Generic;
using System.Linq;

namespace ptmrs
{
    internal class PTMSiteProbability
    {
        internal class SiteProbability
        {
            public ModificationPosition Position;
            public AminoAcidModification Modification;
            public double Probability;

            public SiteProbability(ModificationPosition modPosition, AminoAcidModification modification,
                double probability)
            {
                Position = modPosition;
                Modification = modification;
                Probability = probability;
            }

            public SiteProbability(int position, int numberOfModifications, AminoAcidModification modification,
                double probability)
                : this(new ModificationPosition(position, numberOfModifications, true), modification, probability)
            {
            }

            public override string ToString()
            {
                return $"{Math.Round(Probability * 100.0, 2)}% ({Modification.Abbreviation}, pos: {Position})";
            }

            public override bool Equals(object obj)
            {
                if (obj == null)
                {
                    return false;
                }
                SiteProbability siteProbability = obj as SiteProbability;
                return siteProbability != null && Math.Abs(Probability - siteProbability.Probability) <= 1E-05 &&
                       (siteProbability.Position != null || siteProbability.Position == null) &&
                       (siteProbability.Modification != null || Modification == null) &&
                       ((siteProbability.Position == null && Position == null) ||
                        siteProbability.Position.Equals(Position)) &&
                       ((siteProbability.Modification == null && Modification == null) ||
                        siteProbability.Modification.Equals(siteProbability.Modification));
            }

            public override int GetHashCode()
            {
                return (Position == null ? 0 : Position.GetHashCode()) +
                       (Modification == null ? 0 : Modification.GetHashCode()) + (int) Probability;
            }
        }

        private List<SiteProbability> _probabilityMap;

        public int GetAAModificationCount => _probabilityMap?.Count ?? 0;

        internal PTMSiteProbability(List<SiteProbability> siteProbabilty)
        {
            _probabilityMap = siteProbabilty;
        }

        public override string ToString()
        {
            if (_probabilityMap == null || _probabilityMap.Count <= 0)
            {
                return "";
            }
            return string.Join(";", (
                from w in _probabilityMap
                select w.ToString()).ToArray());
        }

        public string GetSiteProbabilityStrings(string sequence, AminoAcidModification scoredAA)
        {
            return GetSiteProbabilityString(_probabilityMap, sequence, scoredAA);
        }

        private static string GetSiteProbabilityString(List<SiteProbability> probabilityMap,
            string sequence, AminoAcidModification scoredAA)
        {
            string result;
            List<IGrouping<int, SiteProbability>> list = (
                from w in probabilityMap
                where scoredAA.Equals(w.Modification)
                group w by w.Position.ZeroBasedPosition
                into w
                orderby w.Key
                select w).ToList();
            if (list.Count == 0)
                result = "";
            else
            {
                List<string> list2 = new List<string>(list.Sum(w => w.Count()));
                foreach (IGrouping<int, SiteProbability> current in list)
                {
                    int maxNumber = current.Max(w => w.Position.NumberOfModifications);
                    list2.Add(string.Join(";", from s in current
                        orderby s.Position.NumberOfModifications select
                        $"{sequence[s.Position.ZeroBasedPosition]}{s.Position.ToString(maxNumber)}: {Math.Round(s.Probability * 100.0, 2)}"));
                }
                result = string.Join("; ", list2);
            }
            return result;
        }

        public string GetSiteProbabilityString(string sequence,
            List<Tuple<AminoAcidModification, List<AminoAcidSequence>>> aminoAcidModificationIsoformlists)
        {
            string result;
            List<SiteProbability> list =
                aminoAcidModificationIsoformlists.SelectMany(
                        w => from v in _probabilityMap
                            where w.Item2.Any(isoform => isoform.HasModification(v))
                            select v)
                    .Distinct()
                    .ToList();
            List<SiteProbability> source = (
                from w in _probabilityMap.Except(list)
                orderby w.Position.NumberOfModifications
                select w).ToList();
            List<IGrouping<int, SiteProbability>> list2 = (
                from w in list
                orderby w.Position
                group w by w.Position.ZeroBasedPosition).ToList();
            if (list2.Count == 0)
                result = null;
            else
            {
                List<string> list3 =
                    new List<string>(
                        list2.Sum(w =>w.Count()));
                foreach (IGrouping<int, SiteProbability> current in list2)
                {
                    foreach (IGrouping<AminoAcidModification, SiteProbability> current2 in
                        from w in current
                        group w by w.Modification)
                    {
                        int maxnumberofmodifications = current2.Max(w => w.Position.NumberOfModifications);
                        using (IEnumerator<SiteProbability> enumerator3 = current2.GetEnumerator())
                        {
                            while (enumerator3.MoveNext())
                            {
                                SiteProbability prob = enumerator3.Current;
                                var siteProbability = source
                                        .FirstOrDefault(
                                            w => w.Position.OneBasedPosition == prob.Position.OneBasedPosition &&
                                                 w.Position.NumberOfModifications >
                                                 prob.Position.NumberOfModifications &&
                                                 prob.Modification.Equals(w.Modification));
                                double num;
                                if (siteProbability != null)
                                    num = siteProbability.Probability;
                                else
                                    num = 0.0;
                                list3.Add(
                                    $"{sequence[prob.Position.ZeroBasedPosition]}{prob.Position.ToString(maxnumberofmodifications, current2.Key)}: {Math.Round((prob.Probability - num) * 100.0, 2)}");
                            }
                        }
                    }
                }
                result = string.Join("; ", list3);
            }
            return result;
        }
    }
}