using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Xml.Linq;

namespace ptmrs.Data
{
    public class ModificationManagement
    {
        private class EquivalentModificationItem
        {
            public int _NumberOfModifications;
            public char _Target;
            public AminoAcidModification _Modification;

            public double MassDelta
            {
                get
                {
                    if (_Modification == null)
                        return 0.0;
                    return _Modification.MassDelta * _NumberOfModifications;
                }
            }

            public EquivalentModificationItem(int numberOfModifications, AminoAcidModification modification, char target)
            {
                _NumberOfModifications = numberOfModifications;
                _Modification = modification;
                _Target = target;
            }

            public EquivalentModificationItem(Tuple<int, Tuple<AminoAcidModification, char>> modification)
            {
                _NumberOfModifications = modification.Item1;
                _Modification = modification.Item2.Item1;
                _Target = modification.Item2.Item2;
            }

            public override int GetHashCode()
            {
                if (_Modification == null)
                    return base.GetHashCode();
                return _Modification.GetHashCode() + _Target.GetHashCode();
            }

            public override bool Equals(object obj)
            {
                if (_Modification == null)
                    return base.Equals(obj);
                EquivalentModificationItem equivalentModificationItem =
                    obj as EquivalentModificationItem;
                return equivalentModificationItem != null &&
                       equivalentModificationItem.GetHashCode() == GetHashCode() &&
                       equivalentModificationItem._NumberOfModifications == _NumberOfModifications &&
                       equivalentModificationItem._Target.Equals(_Target);
            }

            public override string ToString()
            {
                if (_Modification == null)
                    return base.ToString();
                return $"{_NumberOfModifications}x{_Modification} on {_Target}";
            }
        }

        private Dictionary<Tuple<AminoAcidModification, char>, EquivalentModificationItem> _definedEquivalents;
        private Dictionary<EquivalentModificationItem, Tuple<AminoAcidModification, char>> _definedEquivalents2;

        private List<AminoAcidModification> _scoredModifications;

        public List<AminoAcidModification> ScoredModifications
        {
            get
            {
                if (_scoredModifications == null)
                {
                    if (_definedEquivalents2 == null)
                        return new List<AminoAcidModification>();
                    _scoredModifications = (from equi in _definedEquivalents2
                            select equi.Key._Modification).Distinct()
                        .ToList();
                }
                return _scoredModifications;
            }
        }

        internal Dictionary<AminoAcidModification, List<Tuple<AminoAcid, int>>> m_maxModificationNumber
        {
            get;
            private set;
        }

        public void AddEquivalent(AminoAcidModification modification, char target, int numberOfModifications,
            AminoAcidModification equivalentModi)
        {
            if (_definedEquivalents == null)
                _definedEquivalents = new Dictionary<Tuple<AminoAcidModification, char>, EquivalentModificationItem>();
            if (_definedEquivalents2 == null)
                _definedEquivalents2 = new Dictionary <EquivalentModificationItem, Tuple<AminoAcidModification, char>>();
            EquivalentModificationItem equivalentModificationItem =
                new EquivalentModificationItem(numberOfModifications, equivalentModi, target);
            Tuple<AminoAcidModification, char> tuple = new Tuple<AminoAcidModification, char>(modification, target);

            _definedEquivalents.Add(tuple, equivalentModificationItem);
            _definedEquivalents2.Add(equivalentModificationItem, tuple);
            m_maxModificationNumber = (from w in _definedEquivalents2
                select w.Key
                into w
                group w by w._Modification).ToDictionary(
                v => v.Key,
                v =>(from x in v
                        group x by x._Target
                        into x
                        select (from f in x
                            orderby f._NumberOfModifications descending
                            select f).First()
                        into w
                        select
                        new Tuple<AminoAcid, int>(AminoAcid.GetAminoAcid(w._Target), w._NumberOfModifications))
                    .ToList());
            _scoredModifications = null;
        }

        public Tuple<int, AminoAcidModification> GetEquivalent(AminoAcidModification modification, char target)
        {
            if (_definedEquivalents == null)
                return new Tuple<int, AminoAcidModification>(1, modification);
            Tuple<AminoAcidModification, char> key = new Tuple<AminoAcidModification, char>(modification, target);
            EquivalentModificationItem equivalentModificationItem;
            if (_definedEquivalents.TryGetValue(key, out equivalentModificationItem))
            {
                return new Tuple<int, AminoAcidModification>(equivalentModificationItem._NumberOfModifications,
                    equivalentModificationItem._Modification);
            }
            return new Tuple<int, AminoAcidModification>(1, modification);
        }

        public Tuple<AminoAcidModification, char> GetEquivalent(
            Tuple<int, Tuple<AminoAcidModification, char>> modification)
        {
            if (_definedEquivalents2 == null)
                return modification.Item2;
            EquivalentModificationItem key = new EquivalentModificationItem(modification);
            Tuple<AminoAcidModification, char> result;
            if (_definedEquivalents2.TryGetValue(key, out result))
                return result;
            return modification.Item2;
        }

        public Dictionary<int, AminoAcidModification> GetEquivalent(string sequence,
            Dictionary<int, Tuple<int, AminoAcidModification>> modificationMap)
        {
            if (modificationMap == null)
                return null;
            return modificationMap.ToDictionary(w => w.Key,
                w =>GetEquivalent(new Tuple<int, Tuple<AminoAcidModification, char>>(w.Value.Item1,
                            new Tuple<AminoAcidModification, char>(w.Value.Item2, sequence[w.Key])))
                        .Item1);
        }

        public Dictionary<int, Tuple<int, AminoAcidModification>> GetEquivalent(string sequence,
            Dictionary<int, AminoAcidModification> modificationMap)
        {
            if (modificationMap == null)
                return null;
            return modificationMap.ToDictionary(w => w.Key, w => GetEquivalent(w.Value, sequence[w.Key]));
        }

        public IEnumerable<List<AminoAcidSequence>> GetIsoformPeptideGroups(IEnumerable<AminoAcidSequence> peptides)
        {
            foreach (IGrouping<string, AminoAcidSequence> current in from g in peptides
                group g by g.OneLetterCodeString())
            {
                List<AminoAcidSequence> list =
                    current.OrderBy(o => o.Rank).ToList();
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
                            peptideItem =>
                                peptideItem.GetModificationAASeqPosPairs(false)
                                    .SelectMany(s =>Enumerable.Repeat(s.Value.Item2.ID, s.Value.Item1))
                                    .OrderBy(o => o)
                                    .ToList();
                        Func<AminoAcidSequence, int> getUnusedPTMModifications =
                            peptideItem =>
                                peptideItem
                                    .GetModificationAASeqPosPairs(
                                        false)
                                    .Count(posAndMod => !_definedEquivalents2.Keys.Any(
                                        delegate(EquivalentModificationItem w)
                                        {
                                            if (w._Modification.Equals(posAndMod.Value.Item2) &&
                                                w._NumberOfModifications == posAndMod.Value.Item1)
                                            {
                                                return
                                                    w._Modification.TargetAminoAcids.Select(
                                                            aa => aa.OneLetterCode)
                                                        .Contains(
                                                            peptideItem.GetAminoAcid(posAndMod.Key)
                                                                .OneLetterCode);
                                            }
                                            return false;
                                        }));
                        List<int> modificationsOfCurrentIsoform =
                            getModificationIdsByAscendingOrder(list.First());
                        int numberOfUnusedMods = getUnusedPTMModifications(list.First());
                        List<AminoAcidSequence> list2 =
                            list.FindAll(
                                r =>
                                    getModificationIdsByAscendingOrder(r)
                                        .SequenceEqual(modificationsOfCurrentIsoform) &&
                                    getUnusedPTMModifications(r) == numberOfUnusedMods);
                        list = list.Except(list2).ToList();
                        yield return list2;
                    }
                }
            }
        }

        public AminoAcidSequence GetMappedAminoAcidSequence(AminoAcidSequence aasequence)
        {
            Dictionary<int, Tuple<int, AminoAcidModification>> equivalent =
                this.GetEquivalent(aasequence.OneLetterCodeString(),
                    aasequence.GetModificationAASeqPosPairs(false)
                        .ToDictionary(w => w.Key,
                            w => w.Value.Item2));
            List<AminoAcid> list = (from w in aasequence.OneLetterCodeString()
                select AminoAcid.GetAminoAcid(w)).ToList();
            if (equivalent.Count > 0)
            {
                foreach (KeyValuePair<int, Tuple<int, AminoAcidModification>> current in equivalent)
                {
                    list[current.Key].setModification(current.Value);
                }
            }
            return new AminoAcidSequence(aasequence.ID, aasequence.Rank, list, aasequence.NTerminus,
                aasequence.CTerminus);
        }

        internal void ReadXML_EquivalentModifications(XDocument doc,
            List<AminoAcidModification> possibleAminoAcidModifications)
        {
            foreach (XElement current in doc.Root.Elements())
            {
                if (string.Compare(current.Name.LocalName, "modification", StringComparison.OrdinalIgnoreCase) == 0)
                {
                    string equivalentModificationName = current.Attribute("name").Value;
                    XAttribute xAttribute = current.Attribute("unimodId");
                    AminoAcidModification aminoAcidModification = null;
                    XAttribute xAttribute2 = current.Attribute("searchdefined");
                    if (xAttribute2 == null || string.Compare(xAttribute2.Value, "TRUE", true) == 0)
                    {
                        if (xAttribute != null)
                        {
                            long unimodID = Convert.ToInt64(xAttribute.Value);
                            aminoAcidModification =
                                possibleAminoAcidModifications.Find(w => w.UnimodAccession == unimodID);
                        }
                        if (aminoAcidModification == null)
                        {
                            aminoAcidModification =
                                possibleAminoAcidModifications.Find(
                                    w =>
                                        string.Compare(w.Name, equivalentModificationName, true) == 0 ||
                                        string.Compare(w.Abbreviation, equivalentModificationName, true) == 0);
                        }
                    }
                    else
                    {
                        string value = current.Attribute("abbreviation").Value;
                        double massDelta = Convert.ToDouble(current.Attribute("mass").Value,
                            Constants.EN_US_CULTURE_INFO);
                        List<AminoAcid> list = new List<AminoAcid>();
                        foreach (XElement current2 in
                            current.Elements()
                                .Where(w => string.Compare(w.Name.LocalName, "target", true) == 0))
                        {
                            char letter = Convert.ToChar(current2.Attribute("aminoacid").Value,
                                Constants.EN_US_CULTURE_INFO);
                            list.Add(AminoAcid.GetAminoAcid(letter));
                        }
                        list = list.Distinct().ToList();
                        long unimodAccession = 0L;
                        if (xAttribute != null)
                        {
                            unimodAccession = Convert.ToInt64(xAttribute.Value);
                        }
                        List<DiagnosticIon> diagnosticIons =
                            current.Elements()
                                .Where(w => string.Compare(w.Name.LocalName, "diagnosticion", true) == 0)
                                .Select(w => new DiagnosticIon(w))
                                .ToList();
                        List<Tuple<AminoAcidModification.TargetInformation, AminoAcidModification.NeutralLoss>>
                            list2 =
                                new List
                                <
                                    Tuple
                                    <AminoAcidModification.TargetInformation,
                                        AminoAcidModification.NeutralLoss>>();
                        foreach (XElement current3 in
                            current.Elements()
                                .Where(w => string.Compare(w.Name.LocalName, "neutralloss", true) == 0))
                        {
                            AminoAcidModification.NeutralLoss nl = new AminoAcidModification.NeutralLoss(current3);
                            if (
                                !current3.Elements()
                                    .Where(w => string.Compare(w.Name.LocalName, "target", true) == 0)
                                    .Any())
                            {
                                list2 =
                                    list.Select(
                                            w =>
                                                new Tuple
                                                <AminoAcidModification.TargetInformation,
                                                    AminoAcidModification.NeutralLoss>(
                                                    new AminoAcidModification.TargetInformation(w.OneLetterCode), nl))
                                        .ToList();
                            }
                            else
                            {
                                foreach (XElement current4 in
                                    current3.Elements()
                                        .Where(w => string.Compare(w.Name.LocalName, "target", true) == 0))
                                {
                                    AminoAcid aminoAcid =
                                        AminoAcid.GetAminoAcid(current4.Attribute("aminoacid").Value[0]);
                                    int multiplicity;
                                    if (current4.Attribute("factor") == null)
                                    {
                                        multiplicity = -1;
                                    }
                                    else
                                    {
                                        multiplicity = Convert.ToInt32(current4.Attribute("factor").Value);
                                    }
                                    AminoAcidModification.TargetInformation item =
                                        new AminoAcidModification.TargetInformation(multiplicity,
                                            aminoAcid.OneLetterCode);
                                    list2.Add(
                                        new Tuple
                                        <AminoAcidModification.TargetInformation,
                                            AminoAcidModification.NeutralLoss>(item, nl));
                                }
                            }
                        }
                        if (list2.Count > 0)
                        {
                            Dictionary
                                <AminoAcidModification.TargetInformation, List<AminoAcidModification.NeutralLoss>>
                                neutralLosses = list2.GroupBy(w => w.Item1)
                                    .ToDictionary(w => w.Key, w => w.Select(v => v.Item2).ToList());
                            aminoAcidModification = new AminoAcidModification(AminoAcidModification.getUniqueID(),
                                equivalentModificationName, value, massDelta, list, neutralLosses, diagnosticIons);
                        }
                        else
                        {
                            aminoAcidModification = new AminoAcidModification(AminoAcidModification.getUniqueID(),
                                equivalentModificationName, value, massDelta, list, null, diagnosticIons);
                        }
                        aminoAcidModification.UnimodAccession = unimodAccession;
                    }
                    if (aminoAcidModification != null)
                    {
                        if (
                            current.Elements()
                                .Any(w => string.Compare(w.Name.LocalName, "equivalentmodification", true) == 0))
                        {
                            using (
                                IEnumerator<XElement> enumerator5 =
                                    current.Elements()
                                        .Where(
                                            w =>
                                                string.Compare(w.Name.LocalName, "equivalentmodification", true) ==
                                                0)
                                        .GetEnumerator())
                            {
                                while (enumerator5.MoveNext())
                                {
                                    XElement current5 = enumerator5.Current;
                                    string baseModificationName = current5.Attribute("name").Value;
                                    int num = Convert.ToInt32(current5.Attribute("factor").Value,
                                        Constants.EN_US_CULTURE_INFO);
                                    List<char> list3;
                                    if (current5.Attribute("avoidTarget") != null)
                                    {
                                        list3 =
                                            Convert.ToString(current5.Attribute("avoidTarget").Value,
                                                    Constants.EN_US_CULTURE_INFO)
                                                .ToList<char>();
                                    }
                                    else
                                    {
                                        list3 = null;
                                    }
                                    AminoAcidModification aminoAcidModification2 = null;
                                    XAttribute xAttribute3 = current5.Attribute("unimodId");
                                    if (xAttribute3 != null)
                                    {
                                        long unimodID = Convert.ToInt64(xAttribute3.Value);
                                        aminoAcidModification2 =
                                            possibleAminoAcidModifications.Find(
                                                w => w.UnimodAccession == unimodID);
                                    }
                                    if (aminoAcidModification2 == null)
                                    {
                                        aminoAcidModification2 =
                                            possibleAminoAcidModifications.Find(
                                                w =>
                                                    string.Compare(w.Name, baseModificationName, true) == 0 ||
                                                    string.Compare(w.Abbreviation, baseModificationName, true) == 0);
                                    }
                                    if (aminoAcidModification2 != null)
                                    {
                                        XAttribute xAttribute4 = current5.Attribute("new");
                                        if (xAttribute4 != null &&
                                            string.Compare(xAttribute4.Value, "TRUE", true) == 0)
                                        {
                                            double massDelta2 = Convert.ToDouble(current5.Attribute("mass").Value,
                                                Constants.EN_US_CULTURE_INFO);
                                            string value2 = current5.Attribute("abbreviation").Value;
                                            XElement xElement = current5.Element("neutral");
                                            if (xElement != null)
                                            {
                                                Convert.ToDouble(xElement.Attribute("mass").Value,
                                                    Constants.EN_US_CULTURE_INFO);
                                            }
                                            long unimodAccession2 = 0L;
                                            if (xAttribute3 != null)
                                            {
                                                unimodAccession2 = Convert.ToInt64(xAttribute3.Value);
                                            }
                                            aminoAcidModification2 =
                                                new AminoAcidModification(AminoAcidModification.getUniqueID(),
                                                    baseModificationName, value2, massDelta2, aminoAcidModification);
                                            aminoAcidModification2.UnimodAccession = unimodAccession2;
                                        }
                                        using (
                                            List<AminoAcid>.Enumerator enumerator6 =
                                                aminoAcidModification.TargetAminoAcids.GetEnumerator())
                                        {
                                            while (enumerator6.MoveNext())
                                            {
                                                AminoAcid target = enumerator6.Current;
                                                if (list3 != null)
                                                {
                                                    if (
                                                        list3.Any(
                                                            avoidTarget =>
                                                                avoidTarget.Equals(target.OneLetterCode)))
                                                    {
                                                        continue;
                                                    }
                                                }
                                                AminoAcidModification.TargetInformation targetInformation =
                                                    new AminoAcidModification.TargetInformation(num,
                                                        target.OneLetterCode);
                                                if (aminoAcidModification2.NeutralLossCount(targetInformation) > 0)
                                                {
                                                    foreach (AminoAcidModification.NeutralLoss current6 in
                                                        aminoAcidModification2.NeutralLosses(num,
                                                            target.OneLetterCode))
                                                    {
                                                        aminoAcidModification.AddNeutralLosses(targetInformation,
                                                            current6);
                                                    }
                                                }
                                                AddEquivalent(aminoAcidModification2, target.OneLetterCode, num,
                                                    aminoAcidModification);
                                            }
                                        }
                                    }
                                }
                                continue;
                            }
                        }
                        foreach (AminoAcid current7 in aminoAcidModification.TargetAminoAcids)
                        {
                            new AminoAcidModification.TargetInformation(1, current7.OneLetterCode);
                            AddEquivalent(aminoAcidModification, current7.OneLetterCode, 1,
                                aminoAcidModification);
                        }
                    }
                }
            }
        }

        public IEnumerable<Tuple<char, AminoAcidModification, Tuple<int, AminoAcidModification>>>
            GetScoredEquivalentModifications()
        {
            if (_definedEquivalents2 != null)
            {
                foreach (
                    KeyValuePair<EquivalentModificationItem, Tuple<AminoAcidModification, char>>
                        current in _definedEquivalents2)
                {
                    KeyValuePair<EquivalentModificationItem, Tuple<AminoAcidModification, char>>
                        keyValuePair = current;
                    char argBc0 = keyValuePair.Value.Item2;
                    KeyValuePair<EquivalentModificationItem, Tuple<AminoAcidModification, char>>
                        keyValuePair2 = current;
                    AminoAcidModification argBc1 = keyValuePair2.Value.Item1;
                    KeyValuePair<EquivalentModificationItem, Tuple<AminoAcidModification, char>>
                        keyValuePair3 = current;
                    int argB70 = keyValuePair3.Key._NumberOfModifications;
                    KeyValuePair<EquivalentModificationItem, Tuple<AminoAcidModification, char>>
                        keyValuePair4 = current;
                    yield return
                        new Tuple<char, AminoAcidModification, Tuple<int, AminoAcidModification>>(argBc0, argBc1,
                            new Tuple<int, AminoAcidModification>(argB70, keyValuePair4.Key._Modification));
                }
            }
        }

        internal bool Exists(AminoAcidModification modification, char target)
        {
            bool result;
            if (_definedEquivalents == null)
            {
                result = false;
            }
            else
            {
                Tuple<AminoAcidModification, char> key = new Tuple<AminoAcidModification, char>(modification, target);
                result = _definedEquivalents.ContainsKey(key);
            }
            return result;
        }

        internal bool Exists(AminoAcidModification modification, char target, int numberOfModifications)
        {
            bool result;
            if (_definedEquivalents2 == null)
            {
                result = false;
            }
            else
            {
                EquivalentModificationItem key =
                    new EquivalentModificationItem(numberOfModifications, modification,
                        target);
                result = _definedEquivalents2.ContainsKey(key);
            }
            return result;
        }

        internal bool Exists(AminoAcidModification modification, char target, int numberOfModifications,
            AminoAcidModification equivalentModi)
        {
            bool result;
            Tuple<AminoAcidModification, char> key = new Tuple<AminoAcidModification, char>(equivalentModi, target);
            EquivalentModificationItem obj =
                new EquivalentModificationItem(numberOfModifications, modification, target);
            EquivalentModificationItem equivalentModificationItem;
            result = _definedEquivalents.TryGetValue(key, out equivalentModificationItem) &&
                     equivalentModificationItem.Equals(obj);
            return result;
        }

        internal void CompleteEquivalentModificationAdding(List<AminoAcidModification> scoredAAMods,
            double massTolerance, out List<string> report)
        {
            report = new List<string>();
            if (scoredAAMods != null && scoredAAMods.Count > 0)
            {
                foreach (
                    var current in
                    from aamod in
                    scoredAAMods.SelectMany(w => from aa in w.TargetAminoAcids
                        select new
                        {
                            modification = w,
                            target = aa.OneLetterCode
                        })
                    where !Exists(aamod.modification, aamod.target)
                    select aamod)
                {
                    AddEquivalent(current.modification, current.target, 1, current.modification);
                }
                foreach (
                    IGrouping<char, EquivalentModificationItem> current2 in
                    from w in _definedEquivalents2.Keys
                    group w by w._Target)
                {
                    List<EquivalentModificationItem> list =
                        current2.OrderByDescending(
                                w => w.MassDelta)
                            .ToList();
                    for (int i = 0; i < list.Count; i++)
                    {
                        EquivalentModificationItem equivalentModificationItem = list[i];
                        for (int j = i + 1; j < list.Count; j++)
                        {
                            EquivalentModificationItem equivalentModificationItem2 = list[j];
                            if (
                                Math.Abs(equivalentModificationItem.MassDelta -
                                         equivalentModificationItem2.MassDelta) <= massTolerance)
                            {
                                report.Add(
                                    string.Format("WARNING: Equal amino acid modification detected. {0} ~= {1}.",
                                        equivalentModificationItem._Modification.Name,
                                        equivalentModificationItem2._Modification.Name));
                            }
                        }
                    }
                }
                foreach (
                    Tuple<char, AminoAcidModification, Tuple<int, AminoAcidModification>> current3 in
                    GetScoredEquivalentModifications())
                {
                    AminoAcidModification item = current3.Item3.Item2;
                    int item2 = current3.Item3.Item1;
                    AminoAcidModification item3 = current3.Item2;
                    char item4 = current3.Item1;
                    StringBuilder stringBuilder = new StringBuilder();
                    stringBuilder.Append("Localizing ");
                    stringBuilder.Append(item2);
                    stringBuilder.Append("x");
                    stringBuilder.Append(item.Name);
                    if (!item3.Equals(item))
                    {
                        stringBuilder.Append(" based on ");
                        stringBuilder.Append(item3.Name);
                    }
                    stringBuilder.Append(", elementTarget: ");
                    stringBuilder.Append(item4);
                    stringBuilder.Append(" (delta mass: ");
                    stringBuilder.Append((item.MassDelta * item2).ToString("0.0000"));
                    if (item.NeutralLossCount(item4, item2) == 0)
                    {
                        stringBuilder.Append(", no neutral loss)");
                    }
                    else
                    {
                        stringBuilder.Append(", Neutral-losses:");
                        stringBuilder.Append(item.ToStringNeutralLosses(item4, item2));
                        stringBuilder.Append(")");
                    }
                    if (item.DiagnosticIons != null)
                    {
                        List<string> list2 = new List<string>();
                        foreach (DiagnosticIon current4 in item.DiagnosticIons)
                        {
                            switch (current4.EvidenceType)
                            {
                                case DiagnosticIon.EvidenceTypes.Report:
                                    list2.Add(current4.ToString());
                                    break;
                                case DiagnosticIon.EvidenceTypes.Factor:
                                    if (current4.EvidenceFactor.Contains(item2))
                                    {
                                        list2.Add(current4.ToString());
                                    }
                                    break;
                                case DiagnosticIon.EvidenceTypes.Target:
                                    if (current4.EvidenceAminoAcid.Contains(item4))
                                    {
                                        list2.Add(current4.ToString());
                                    }
                                    break;
                                case DiagnosticIon.EvidenceTypes.TargetAndFactor:
                                    if (current4.EvidenceAminoAcid.Contains(item4) &&
                                        current4.EvidenceFactor.Contains(item2))
                                    {
                                        list2.Add(current4.ToString());
                                    }
                                    break;
                            }
                        }
                        if (list2.Count > 0)
                        {
                            stringBuilder.Append(", diagnostic-ion(s): ");
                            stringBuilder.Append(string.Join(", ", list2));
                        }
                    }
                    report.Add(stringBuilder.ToString());
                }
            }
        }
    }
}