using System;

namespace ptmrs
{
    // 2017-04-19
    public class FragmentIon : IComparable<double>
    {
        public double Mass;
        public FragmentIonType Type;
        public int Index;

        /// <summary>
        /// Constructor.
        /// </summary>
        /// <param name="mass">离子质量</param>
        /// <param name="type">离子类型</param>
        /// <param name="index">ion index</param>
        public FragmentIon(double mass, FragmentIonType type, int index)
        {
            Mass = mass;
            Type = type;
            Index = index;
        }

        public int CompareTo(double other)
        {
            return other.CompareTo(Mass);
        }

        public override string ToString()
        {
            return $"{Type}({Index}): {Mass:F2}";
        }

        public override bool Equals(object obj)
        {
            if (obj == null)
                return false;
            FragmentIon fragmentIon = obj as FragmentIon;
            return fragmentIon != null && PtmMathHelper.Equal(Mass, fragmentIon.Mass, 0.0001) &&
                   Type.Equals(fragmentIon.Type) && Index == fragmentIon.Index;
        }

        public override int GetHashCode()
        {
            return (int) Mass + Type.GetHashCode() + Index;
        }
    }
}