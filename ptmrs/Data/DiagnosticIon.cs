using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Xml.Linq;

namespace ptmrs.Data
{
    public class DiagnosticIon
    {
        public enum EvidenceTypes
        {
            Report,
            Factor,
            Target,
            TargetAndFactor
        }

        [Flags]
        public enum Flaqs
        {
            NotSet = 0,
            SingleChargedOnly = 1,
            AllChargeStates = 2
        }

        public enum PeakType
        {
            ExistingPrecursor,
            Precursor,
            ImoniumIon
        }

        private int m_peakdepth;


        private readonly PeakType _relativity;

        private readonly SpectrumType _spectrumType;

        private DiagnosticIon(string name, double mass, int peakDepth, SpectrumType spectrumType)
        {
            this.Name = name;
            this.Mass = mass;
            this.PeakDepth = peakDepth;
            EvidenceFactor = new List<int>(1) { 1 };
            EvidenceType = EvidenceTypes.Report;
            _relativity = PeakType.ExistingPrecursor;
            _spectrumType = spectrumType;
        }

        public DiagnosticIon(string name, double mass, int peakDepth, List<int> factor, SpectrumType spectrumType)
            : this(name, mass, peakDepth, spectrumType)
        {
            EvidenceFactor = factor;
            EvidenceType = EvidenceTypes.Report;
        }

        public DiagnosticIon(string name, double mass, int peakDepth, List<char> target, SpectrumType spectrumType)
            : this(name, mass, peakDepth, spectrumType)
        {
            EvidenceAminoAcid = target;
            EvidenceType = EvidenceTypes.Target;
        }

        public DiagnosticIon(string name, double mass, int peakDepth, List<char> target, List<int> factor,
            SpectrumType spectrumType) : this(name, mass, peakDepth, spectrumType)
        {
            EvidenceFactor = factor;
            EvidenceAminoAcid = target;
            EvidenceType = EvidenceTypes.TargetAndFactor;
        }

        public DiagnosticIon(XElement diagnosticXML)
            : this(
                diagnosticXML.Attribute("name").Value,
                Convert.ToDouble(diagnosticXML.Attribute("mass").Value, Constants.EN_US_CULTURE_INFO),
                Convert.ToInt32(diagnosticXML.Attribute("peakdepth").Value, Constants.EN_US_CULTURE_INFO),
                SpectrumType.None)
        {
            if (diagnosticXML.Element("spectrumType") == null)
            {
                _spectrumType = SpectrumType.None;
            }
            else if (
                !Enum.TryParse(diagnosticXML.Element("spectrumType").Value, true,
                    out _spectrumType))
            {
                _spectrumType = SpectrumType.None;
            }
            if (diagnosticXML.Element("Relativity") == null)
            {
                _relativity = PeakType.ExistingPrecursor;
            }
            else
            {
                var value = diagnosticXML.Element("Relativity").Value;
                if (string.Compare(value, "ExistingPrecursor", true) == 0)
                {
                    _relativity = PeakType.ExistingPrecursor;
                }
                else if (string.Compare(value, "Precursor", true) == 0)
                {
                    _relativity = PeakType.Precursor;
                }
                else
                {
                    _relativity = PeakType.ImoniumIon;
                }
            }
            var flag = false;
            var flag2 = false;
            if (diagnosticXML.Elements("Evidence_Factor").Count() > 0)
            {
                flag2 = true;
                EvidenceFactor = new List<int>();
                foreach (var current in diagnosticXML.Elements("Evidence_Factor"))
                {
                    EvidenceFactor.Add(Convert.ToInt32(current.Value, Constants.EN_US_CULTURE_INFO));
                }
            }
            if (diagnosticXML.Elements("Evidence_Target").Count() > 0)
            {
                flag = true;
                EvidenceAminoAcid = new List<char>();
                foreach (var current2 in diagnosticXML.Elements("Evidence_Target"))
                {
                    EvidenceAminoAcid.Add(Convert.ToChar(current2.Value, Constants.EN_US_CULTURE_INFO));
                }
            }
            if (flag && flag2)
            {
                EvidenceType = EvidenceTypes.TargetAndFactor;
                return;
            }
            if (flag)
            {
                EvidenceType = EvidenceTypes.Target;
                return;
            }
            if (flag2)
            {
                EvidenceType = EvidenceTypes.Factor;
                return;
            }
            throw new ArgumentException("Evidence type not set");
        }

        public double Mass { get; }

        public string Name { get; }

        public int PeakDepth { get; }

        public EvidenceTypes EvidenceType { get; }

        public List<int> EvidenceFactor { get; }

        public List<char> EvidenceAminoAcid { get; }

        public SpectrumType SpectrumType => _spectrumType;

        public override string ToString()
        {
            string result;
            try
            {
                var stringBuilder = new StringBuilder();
                stringBuilder.Append(Name);
                stringBuilder.Append(" (mass: ");
                stringBuilder.Append(Mass.ToString("0.0000"));
                stringBuilder.Append(", peak depth: ");
                stringBuilder.Append(PeakDepth);
                if (_spectrumType != SpectrumType.None)
                {
                    stringBuilder.Append(", spectrum type: ");
                    stringBuilder.Append(Enum.GetName(typeof(SpectrumType), _spectrumType));
                }
                switch (EvidenceType)
                {
                    case EvidenceTypes.Report:
                        stringBuilder.Append(", report)");
                        break;
                    case EvidenceTypes.Factor:
                        stringBuilder.Append(", evidence factor: ");
                        stringBuilder.Append(string.Join(",", EvidenceFactor));
                        stringBuilder.Append(")");
                        break;
                    case EvidenceTypes.Target:
                        stringBuilder.Append(", evidence target: ");
                        stringBuilder.Append(string.Join(",", EvidenceAminoAcid));
                        stringBuilder.Append(")");
                        break;
                    case EvidenceTypes.TargetAndFactor:
                        stringBuilder.Append(", evidence factor: ");
                        stringBuilder.Append(string.Join(",", EvidenceFactor));
                        stringBuilder.Append(", evidence target: ");
                        stringBuilder.Append(string.Join(",", EvidenceAminoAcid));
                        stringBuilder.Append(")");
                        break;
                }
                result = stringBuilder.ToString();
            }
            catch (Exception)
            {
                result = base.ToString();
            }
            return result;
        }

        private IEnumerable<Tuple<int, double>> GetBasePeak(double precursorMZ, int charge, Flaqs flaqs)
        {
            if (_relativity == PeakType.ImoniumIon)
            {
                yield return new Tuple<int, double>(1, 0.0);
            }
            else if (flaqs == Flaqs.AllChargeStates)
            {
                for (var i = charge; i >= 1; i--)
                {
                    yield return
                        new Tuple<int, double>(i, (precursorMZ * charge + 1.007276452 * (i - charge)) / i);
                }
            }
            else if (flaqs == Flaqs.SingleChargedOnly)
            {
                yield return new Tuple<int, double>(1, precursorMZ * charge - (1 - charge) * 1.007276452);
            }
            else
            {
                yield return new Tuple<int, double>(charge, precursorMZ);
            }
        }

        internal List<AminoAcidSequence> FilterPotentialIsoforms(List<AminoAcidSequence> isoforms,
            int numberOfAAModOnSequence, AminoAcidModification parentModification, int charge, double masstolerance,
            double precursorMZ, PeakExtractor peakExtractor, out string report)
        {
            List<AminoAcidSequence> result;
            report = null;
            var mostIntensePeaks = peakExtractor.GetMostIntensePeaks(PeakDepth);
            int num;
            if (_relativity != PeakType.ImoniumIon)
            {
                num = (int)Math.Floor(numberOfAAModOnSequence / (double)EvidenceFactor.Min());
            }
            else
            {
                num = 1;
            }
            var nMax = 0;
            var list = new List<AminoAcidSequence>(isoforms);
            using (var enumerator = GetBasePeak(precursorMZ, charge, Flaqs.SingleChargedOnly).GetEnumerator())
            {
                double mass;
                while (enumerator.MoveNext())
                {
                    var basePeak = enumerator.Current;
                    var i = 1;
                    while (i <= num)
                    {
                        if (_relativity != PeakType.ExistingPrecursor)
                        {
                            goto IL_F0;
                        }
                        var peak =
                            mostIntensePeaks.Find(p => Math.Abs(basePeak.Item2 - p.MassZ) <= masstolerance);
                        if (peak != null)
                        {
                            goto IL_F0;
                        }
                        IL_150:
                        i++;
                        continue;
                        IL_F0:
                        mass = basePeak.Item2 + i * Mass / basePeak.Item1;
                        peak = mostIntensePeaks.Find(p => Math.Abs(mass - p.MassZ) <= masstolerance);
                        if (peak != null)
                        {
                            nMax = Math.Max(i, nMax);
                        }
                        goto IL_150;
                    }
                    if (nMax == num)
                    {
                        break;
                    }
                }
            }
            if (nMax <= 0)
            {
                result = isoforms;
            }
            else
            {
                switch (EvidenceType)
                {
                    case EvidenceTypes.Report:
                        report = $"{nMax}x{Name}:Report";
                        break;
                    case EvidenceTypes.Factor:
                        report = $"{nMax}x{Name}:Filter, factor {string.Join(",", EvidenceFactor)}";
                        list = (from aaseq in list
                                where
                                    aaseq.GetModificationAASeqPosPairs(false)
                                        .Count(pos => pos.Value.Item2.Equals(parentModification) &&
                                                EvidenceFactor.Contains(pos.Value.Item1)) >= nMax
                                select aaseq).ToList();
                        break;
                    case EvidenceTypes.Target:
                        report = $"{nMax}x{Name}:Filter, target {string.Join(",", EvidenceAminoAcid)}";
                        list = (from aaseq in list
                                where aaseq.GetModificationAASeqPosPairs(false)
                                        .Count(pos =>
                                                pos.Value.Item2.Equals(parentModification) &&
                                                EvidenceAminoAcid.Contains(
                                                    aaseq.GetAminoAcid(pos.Key).OneLetterCode)) >= nMax
                                select aaseq).ToList();
                        break;
                    case EvidenceTypes.TargetAndFactor:
                        report =
                            $"{nMax}x{Name}:Filter, factor {string.Join(",", EvidenceFactor)} and target {string.Join(",", EvidenceAminoAcid)}";
                        list = (from aaseq in list
                                where
                                    aaseq.GetModificationAASeqPosPairs(false)
                                        .Count(pos =>
                                                pos.Value.Item2.Equals(parentModification) &&
                                                EvidenceFactor.Contains(pos.Value.Item1) &&
                                                EvidenceAminoAcid.Contains(aaseq.GetAminoAcid(pos.Key).OneLetterCode)) >= nMax
                                select aaseq).ToList();
                        break;
                }
                result = list;
            }
            return result;
        }
    }
}