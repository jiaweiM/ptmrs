using System.Collections.Generic;

namespace ptmrs
{
    // 2017-4-20
   /// <summary>
   /// 指定离子类型及其所有的离子质量
   /// </summary>
    public class FitMasses
    {
        public const double NoMass = -1.0;

        public FragmentIonType Type;
        public List<FragmentIon> FragmentIons;

        public FitMasses()
        {
            Type = new FragmentIonType();
            FragmentIons = new List<FragmentIon>();
        }

        public FitMasses(FragmentIonType type)
        {
            Type = type;
            FragmentIons = new List<FragmentIon>();
        }

        public FitMasses(FragmentIonType type, List<FragmentIon> fragmentIons)
        {
            Type = type;
            FragmentIons = fragmentIons;
        }
    }
}