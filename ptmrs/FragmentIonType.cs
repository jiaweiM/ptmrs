using System;
using System.Collections.Generic;
using System.Linq;

namespace ptmrs
{
    public class FragmentIonType
    {
        public const int FullExpansion = -1;
        protected FragmentIonTerminalType Terminal;
        protected bool IsMultiNl;
        protected int _expansion;
        protected int _charge;
        protected double _massDiff;
        protected string _title;
        protected string _appendix;

        internal FitComposition FitComposition { get; }

        public FragmentIonTerminalType TerminalType
        {
            get => Terminal;
            set => Terminal = value;
        }

        public int Expansion
        {
            get => _expansion;
            set => _expansion = value;
        }

        public int Charge
        {
            get => _charge;
            set => _charge = value;
        }

        public double MassDiff
        {
            get => _massDiff;
            set => _massDiff = value;
        }

        public string Title
        {
            get => _title;
            set => _title = value;
        }

        public string Appendix
        {
            get => _appendix;
            set => _appendix = value;
        }

        public bool IsNTerminal => Terminal == FragmentIonTerminalType.NTerminal;

        public bool IsCTerminal => Terminal == FragmentIonTerminalType.CTerminal;

        public bool IsNonTerminal => Terminal == FragmentIonTerminalType.NonTerminal;

        public FragmentIonType()
        {
            Terminal = FragmentIonTerminalType.NonTerminal;
            IsMultiNl = false;
            _expansion = -1;
            _charge = 0;
            _massDiff = 0.0;
            _title = "";
            _appendix = "";
        }

        public FragmentIonType(FragmentIonTerminalType terminal, bool isMultiNL, int expansion, int charge,
            double massDiff, string title, string appendices)
        {
            Terminal = terminal;
            IsMultiNl = isMultiNL;
            _expansion = expansion;
            _charge = charge;
            _massDiff = massDiff;
            _title = title;
            _appendix = appendices;
        }

        /// <summary>
        /// 
        /// </summary>
        /// <param name="terminal"></param>
        /// <param name="isMultiNl"></param>
        /// <param name="expansion"></param>
        /// <param name="charge"></param>
        /// <param name="massDiff"></param>
        /// <param name="title"></param>
        /// <param name="appendices"></param>
        public FragmentIonType(FragmentIonTerminalType terminal, bool isMultiNl, int expansion, int charge,
            double massDiff, FitComposition title, string appendices)
        {
            Terminal = terminal;
            IsMultiNl = isMultiNl;
            _expansion = expansion;
            _charge = charge;
            _massDiff = massDiff;
            FitComposition = title;
            if (title == FitComposition.ZPrime || title == FitComposition.ZRadical)
                _title = "z";
            else
                _title = Enum.GetName(typeof(FitComposition), title);
            _appendix = appendices;
        }

        public FragmentIonType(FragmentIonType baseType, int expansion, double massDiff, string appendix,
            string sep = "-")
            : this(baseType.TerminalType, baseType.IsMultiNeutralLossType(), expansion, baseType.Charge,
                baseType.MassDiff + massDiff, baseType.FitComposition, string.Join(sep, new string[]
                    {baseType.Appendix, appendix}))
        {
        }

        public FragmentIonType(FragmentIonType baseType, bool isMultiNL, int expansion, double massDiff,
            string appendix, string sep = "-") : this(baseType, expansion, massDiff, appendix, sep)
        {
            IsMultiNl = isMultiNL;
        }

        public static List<FragmentIonType> GetAllBaseTypes(List<FragmentIonType> allFITs)
        {
            return (from fit in allFITs where fit.IsBaseType() select fit).ToList();
        }

        public static FragmentIonType[] GetAllNLTypes(bool ignorePriority, bool ignoreExpansion, bool ignoreAppendix,
            List<FragmentIonType> allFITs)
        {
            List<FragmentIonType> list = new List<FragmentIonType>();
            foreach (FragmentIonType current in allFITs)
            {
                if (current.IsNeutralLossType())
                {
                    int num = current.Contains(ignorePriority, ignoreExpansion, ignoreAppendix, list);
                    if (num == -1)
                    {
                        list.Add(current);
                    }
                    else
                    {
                        FragmentIonType fragmentIonType = list[num];
                        if ((current.IsNTerminal && current.Expansion < fragmentIonType.Expansion) ||
                            (current.IsCTerminal && current.Expansion > fragmentIonType.Expansion))
                        {
                            list[num] = current;
                        }
                    }
                }
            }
            return list.ToArray();
        }

        public bool IsZIon()
        {
            return _title.Equals("z");
        }

        public bool IsBaseTypeZIon()
        {
            return IsZIon() && (string.IsNullOrEmpty(_appendix) || _appendix.Equals("+H") ||
                    _appendix.Equals("+2H"));
        }

        public bool IsBaseType()
        {
            return !IsNonTerminal && (string.IsNullOrEmpty(Appendix) || IsBaseTypeZIon());
        }

        public bool IsNeutralLossType()
        {
            return !IsNonTerminal && !string.IsNullOrEmpty(Appendix) && !IsBaseTypeZIon();
        }

        public bool IsMultiNeutralLossType()
        {
            return IsMultiNl;
        }

        public bool IsPrecursorType()
        {
            return IsNonTerminal && Expansion == -1;
        }

        public override bool Equals(object obj)
        {
            return Equals(obj, true, false, true, false, true);
        }

        public bool Equals(object obj, bool ignorePriority, bool ignoreExpansion, bool ignoreAppendix, bool ignoreMassDiff, bool ignoreColor)
        {
            FragmentIonType fragmentIonType = obj as FragmentIonType;
            return fragmentIonType != null && (Terminal == fragmentIonType.TerminalType &&
                                               IsMultiNl == fragmentIonType.IsMultiNeutralLossType() &&
                                               (ignoreExpansion || _expansion == fragmentIonType.Expansion) &&
                                               _charge == fragmentIonType.Charge) && (ignoreMassDiff || PtmMathHelper.Equal(_massDiff, fragmentIonType.MassDiff, 1E-06)) && _title.Equals(fragmentIonType.Title) && (ignoreAppendix || _appendix.Equals(fragmentIonType.Appendix));
        }

        public override int GetHashCode()
        {
            double num = Terminal.GetHashCode();
            num += IsMultiNl.GetHashCode();
            num += _charge.GetHashCode();
            num += _title.GetHashCode();
            return (int) (num % int.MaxValue);
        }

        public override string ToString()
        {
            return string.Join("", Title, Charge, "+", Appendix ?? "");
        }
    }
}