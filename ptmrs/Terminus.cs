using System;
using System.Text;

namespace ptmrs
{
    // 4.19
    public class Terminus : AAModTarget
    {
        private FragmentIonTerminalType _type;
        private string _toString;

        public FragmentIonTerminalType Type
        {
            get => _type;
            set => _type = value;
        }

        public Terminus() : this(FragmentIonTerminalType.NTerminal, null)
        {
        }

        public Terminus(FragmentIonTerminalType type) : this(type, null)
        {
        }

        public Terminus(FragmentIonTerminalType type, AminoAcidModification mod)
        {
            _type = type;
            m_mod = mod;
            m_numberOfModifications = 1;
        }

        public override void setModification(AminoAcidModification modification, int number)
        {
            m_mod = modification;
            m_numberOfModifications = number;
        }

        public override void setModification(Tuple<int, AminoAcidModification> modification)
        {
            setModification(modification.Item2, modification.Item1);
        }

        // 4.19
        public override double GetTotalMass()
        {
            double num = 0.0;
            switch (_type)
            {
                case FragmentIonTerminalType.NTerminal:
                    num += Atoms.HydrogenMass;
                    break;
                case FragmentIonTerminalType.CTerminal:
                    num += Atoms.OxygenMass + Atoms.HydrogenMass;
                    break;
            }
            if (m_mod != null)
            {
                num += m_mod.MassDelta;
            }
            return num;
        }

        public override bool Equals(object obj)
        {
            Terminus terminus = (Terminus) obj;
            return terminus != null && _type == terminus.Type &&
                   (m_mod == null && terminus.Modification == null ||
                    m_mod != null && m_mod.Equals(terminus.Modification));
        }

        public override string ToString()
        {
            if (_toString == null)
            {
                StringBuilder stringBuilder = new StringBuilder();
                if (_type == FragmentIonTerminalType.NonTerminal)
                    stringBuilder.Append("Non-Term");
                else if (_type == FragmentIonTerminalType.CTerminal)
                    stringBuilder.Append("C");
                else
                    stringBuilder.Append("N");
                if (IsModified)
                {
                    stringBuilder.Append(" (");
                    stringBuilder.Append(m_mod);
                    stringBuilder.Append("x");
                    stringBuilder.Append(m_numberOfModifications);
                    stringBuilder.Append(")");
                }
                _toString = stringBuilder.ToString();
            }
            return _toString;
        }

        public override int GetHashCode()
        {
            double num = _type.GetHashCode();
            if (m_mod != null)
            {
                num += m_mod.GetHashCode();
            }
            return (int) num % int.MaxValue;
        }
    }
}