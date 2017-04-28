using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace ptmrs
{
    // 4.19
    public class AminoAcid : AAModTarget
    {
        protected char Letter;
        protected string ShortName;
        protected string MName;
        protected double MonoisoMass;
        public static Dictionary<char, AminoAcid> AllAaRes;

        private double _mAaMass = -9999.0;
        private int _myhashCode;

        public bool HasNL => IsModified && (m_mod.AllowsNL(OneLetterCode, m_numberOfModifications) ||
                                            m_mod.AllowsNL(OneLetterCode, -1));

        public char OneLetterCode => Letter;
        public string ThreeLetterCode => ShortName;
        public string Name => MName;
        public double MonoisotopicMass => MonoisoMass;

        static AminoAcid()
        {
            AllAaRes = new Dictionary<char, AminoAcid>
            {
                {'A', new AminoAcid('A', "Ala", "Alanine", 71.037114, null, 0)},
                {'B', new AminoAcid('B', "Asx", "Asn or Asp", 114.53494, null, 0)},
                {'C', new AminoAcid('C', "Cys", "Cysteine", 103.009185, null, 0)},
                {'D', new AminoAcid('D', "Asp", "Aspartic acid", 115.026943, null, 0)},
                {'E', new AminoAcid('E', "Glu", "Glutamic acid", 129.042593, null, 0)},
                {'F', new AminoAcid('F', "Phe", "Phenylalanine", 147.068414, null, 0)},
                {'G', new AminoAcid('G', "Gly", "Glycine", 57.021464, null, 0)},
                {'H', new AminoAcid('H', "His", "Histidine", 137.058912, null, 0)},
                {'I', new AminoAcid('I', "Ile", "Isoleucine", 113.084064, null, 0)},
                {'K', new AminoAcid('K', "Lys", "Lysine", 128.094963, null, 0)},
                {'L', new AminoAcid('L', "Leu", "Leucine", 113.084064, null, 0)},
                {'M', new AminoAcid('M', "Met", "Methionine", 131.040485, null, 0)},
                {'N', new AminoAcid('N', "Asn", "Asparagine", 114.042927, null, 0)},
                {'P', new AminoAcid('P', "Pro", "Proline", 97.052764, null, 0)},
                {'O', new AminoAcid('O', "Pyl", "Pyrrolysine", 237.14772, null, 0)},
                {'Q', new AminoAcid('Q', "Gln", "Glutamine", 128.058578, null, 0)},
                {'R', new AminoAcid('R', "Arg", "Arginine", 156.101111, null, 0)},
                {'S', new AminoAcid('S', "Ser", "Serine", 87.032028, null, 0)},
                {'T', new AminoAcid('T', "Thr", "Threonine", 101.047679, null, 0)},
                {'U', new AminoAcid('U', "SeC", "Selenocysteine", 150.95363, null, 0)},
                {'V', new AminoAcid('V', "Val", "Valine", 99.068414, null, 0)},
                {'W', new AminoAcid('W', "Trp", "Tryptophan", 186.079313, null, 0)},
                {'X', new AminoAcid('X', "Xaa", "Unknown", 111.0, null, 0)},
                {'Y', new AminoAcid('Y', "Tyr", "Tyrosine", 163.06333, null, 0)},
                {'Z', new AminoAcid('Z', "Glx", "Glu or Gln", 128.55059, null, 0)}
            };
        }

        public AminoAcid(char letter, string shortName, string name, double monoisotopicMass,
            AminoAcidModification modification, int numberOfModifications)
        {
            Letter = letter;
            ShortName = shortName;
            MName = name;
            MonoisoMass = monoisotopicMass;
            m_mod = modification;
            m_numberOfModifications = numberOfModifications;
            double num = Letter.GetHashCode();
            _myhashCode = (int) num % int.MaxValue;
        }

        public override void setModification(AminoAcidModification modification, int number)
        {
            m_mod = modification;
            m_numberOfModifications = number;
            _myhashCode = Letter.GetHashCode() % int.MaxValue;
            _mAaMass = -9999.0;
        }

        public override void setModification(Tuple<int, AminoAcidModification> modification)
        {
            setModification(modification.Item2, modification.Item1);
        }

        public static bool FindAAResidue(char letter, out AminoAcid outResult)
        {
            AminoAcid aminoAcid;
            if (AllAaRes.TryGetValue(letter, out aminoAcid))
            {
                outResult = new AminoAcid(aminoAcid.OneLetterCode, aminoAcid.ThreeLetterCode, aminoAcid.Name,
                    aminoAcid.MonoisotopicMass, aminoAcid.Modification, aminoAcid.NumberOfModifications);
                return true;
            }
            outResult = null;
            return false;
        }

        public static AminoAcid GetAminoAcid(char letter)
        {
            AminoAcid result;
            if (!FindAAResidue(letter, out result))
                throw new ArgumentException("Amino Acid not found!");
            return result;
        }

        public static bool FindAAResidue(char letter, AminoAcidModification modification, int numberOfModifications,
            out AminoAcid outResult)
        {
            AminoAcid aminoAcid;
            if (AllAaRes.TryGetValue(letter, out aminoAcid))
            {
                outResult = new AminoAcid(aminoAcid.OneLetterCode, aminoAcid.ThreeLetterCode, aminoAcid.Name,
                    aminoAcid.MonoisotopicMass, modification, numberOfModifications);
                return true;
            }
            outResult = null;
            return false;
        }

        public static double SummAAMasses(List<AminoAcid> aas)
        {
            return SumAAMasses(aas.Count, aas, false);
        }

        public static double SumAAMasses(int cnt, List<AminoAcid> aas, bool reverse)
        {
            double num = 0.0;
            if (cnt <= 0)
                return 0.0;
            if (cnt > aas.Count)
                cnt = aas.Count;
            if (reverse)
            {
                for (int i = aas.Count - cnt; i < aas.Count; i++)
                    num += aas[i].GetTotalMass();
            }
            else
            {
                for (int j = 0; j < cnt; j++)
                    num += aas[j].GetTotalMass();
            }
            return num;
        }

        public static string ParseUnmodified1LetterCodeString(List<AminoAcid> aas)
        {
            StringBuilder stringBuilder = new StringBuilder();
            for (int i = 0; i < aas.Count; i++)
                stringBuilder.Append(aas[i].OneLetterCode);
            return stringBuilder.ToString();
        }

        public override string ToString()
        {
            return (IsModified ? char.ToLower(OneLetterCode) : char.ToUpper(OneLetterCode)).ToString();
        }

        public string GetModificationString()
        {
            if (!IsModified)
                return "";
            if (NumberOfModifications > 1)
                return $"({NumberOfModifications}x{m_mod.Abbreviation.First()})";
            return $"({m_mod.Abbreviation.First()})";
        }

        public override double GetTotalMass()
        {
            if (_mAaMass == -9999.0)
            {
                _mAaMass = MonoisoMass + (IsModified ? m_mod.MassDelta * NumberOfModifications : 0.0);
            }
            return _mAaMass;
        }

        public override bool Equals(object obj)
        {
            return Equals(obj, true);
        }

        public bool Equals(object obj, bool ignoreMod)
        {
            AminoAcid aminoAcid = obj as AminoAcid;
            if (aminoAcid == null || GetHashCode() != aminoAcid.GetHashCode())
                return false;
            if (ignoreMod || m_mod == null && aminoAcid.Modification == null)
                return true;
            if (m_mod != null && aminoAcid.NumberOfModifications == NumberOfModifications)
                return m_mod.Equals(aminoAcid.Modification);
            return false;
        }

        public override int GetHashCode()
        {
            return _myhashCode;
        }
    }
}