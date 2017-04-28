using System;
using System.Collections.Generic;
using System.Linq;

namespace ptmrs
{
    internal class ModificationPosition : IComparable, IComparable<ModificationPosition>
    {
        internal class Compare : IEqualityComparer<List<ModificationPosition>>
        {
            public bool Equals(List<ModificationPosition> x, List<ModificationPosition> y)
            {
                return List_Comparer(x, y);
            }

            public int GetHashCode(List<ModificationPosition> codeh)
            {
                if (codeh == null)
                {
                    return 0;
                }
                return codeh.GetHashCode();
            }
        }

        private readonly int _position;
        public readonly int NumberOfModifications;
        public readonly bool BasedOnZero;

        public int ZeroBasedPosition
        {
            get
            {
                if (!BasedOnZero)
                    return _position - 1;
                return _position;
            }
        }

        public int OneBasedPosition
        {
            get
            {
                if (!BasedOnZero)
                    return _position;
                return _position + 1;
            }
        }

        public ModificationPosition(int position, int numberOfModifications, bool basedOnZero)
        {
            _position = position;
            NumberOfModifications = numberOfModifications;
            BasedOnZero = basedOnZero;
        }

        public override string ToString()
        {
            return ToString(NumberOfModifications);
        }

        public string ToString(int maxnumberofmodifications)
        {
            return ToString(maxnumberofmodifications, null);
        }

        public string ToString(int maxnumberofmodifications, AminoAcidModification modification)
        {
            if (modification == null)
            {
                return $"({OneBasedPosition})x{NumberOfModifications}";
            }
            if (maxnumberofmodifications <= 1)
            {
                return $"{OneBasedPosition}({modification.Abbreviation})";
            }
            if (maxnumberofmodifications == NumberOfModifications)
            {
                return $"{OneBasedPosition}({NumberOfModifications}x{modification.Abbreviation})";
            }
            return $"{OneBasedPosition}(min {NumberOfModifications}x{modification.Abbreviation})";
        }

        public override int GetHashCode()
        {
            return 1000 * _position.GetHashCode() + NumberOfModifications.GetHashCode();
        }

        public override bool Equals(object obj)
        {
            if (obj == null)
            {
                return false;
            }
            ModificationPosition modificationPosition = obj as ModificationPosition;
            return modificationPosition != null &&
                   (modificationPosition._position == _position &&
                    modificationPosition.NumberOfModifications == NumberOfModifications) &&
                   modificationPosition.BasedOnZero == BasedOnZero;
        }

        public static bool List_Comparer(List<ModificationPosition> value1, List<ModificationPosition> value2)
        {
            bool result;
            try
            {
                if (value1 == value2)
                {
                    result = true;
                }
                else if (value1.Count == value2.Count)
                {
                    result = true;
                }
                else
                {
                    result = (!value1.Intersect(value2).Any());
                }
            }
            catch (Exception)
            {
                result = false;
            }
            return result;
        }

        public int CompareTo(object obj)
        {
            int result;
            if (obj == null)
            {
                result = -1;
            }
            else
            {
                ModificationPosition modificationPosition = obj as ModificationPosition;
                if (modificationPosition == null)
                {
                    throw new ArgumentException();
                }
                result = CompareTo(modificationPosition);
            }
            return result;
        }

        public int CompareTo(ModificationPosition other)
        {
            if (other == null)
            {
                return -1;
            }
            if (other.ZeroBasedPosition < ZeroBasedPosition)
            {
                return 1;
            }
            if (other.ZeroBasedPosition > ZeroBasedPosition)
            {
                return -1;
            }
            return NumberOfModifications - other.NumberOfModifications;
        }
    }
}