using System;
using System.Collections.Generic;

// 2017-04-19
namespace ptmrs
{
    public class Peak : IComparable<Peak>
    {
        private Dictionary<AminoAcidSequence, int> _matchingAaSequences;

        public double MassZ { get; set; }

        public double Intensity { get; set; }

        public Peak() : this(0.0, 0.0)
        {
        }

        public Peak(double massZ, double intensity)
        {
            MassZ = massZ;
            Intensity = intensity;
        }

        /// <summary>
        /// Parse a string into  peaks.
        /// </summary>
        /// <param name="s">peaks string representation</param>
        /// <returns>list of Peak</returns>
        public static List<Peak> ParsePeaks(string s)
        {
            string[] array = s.Split(':', ',');
            List<Peak> list = new List<Peak>(array.Length / 2);
            for (int i = 0; i < array.Length / 2; i++)
            {
                double massZ;
                double intensity;
                if (double.TryParse(array[i * 2], out massZ) && double.TryParse(array[i * 2 + 1], out intensity))
                {
                    list.Add(new Peak(massZ, intensity));
                }
            }
            return list;
        }

        public int IsMatchingAASequence(AminoAcidSequence seq)
        {
            int num;
            if (_matchingAaSequences == null || !_matchingAaSequences.TryGetValue(seq, out num))
                return 0;
            return num;
        }

        public void addMatchingAASequence(AminoAcidSequence match)
        {
            if (_matchingAaSequences == null)
            {
                _matchingAaSequences = new Dictionary<AminoAcidSequence, int>();
            }
            int num;
            if (_matchingAaSequences.TryGetValue(match, out num))
                _matchingAaSequences[match] = num + 1;
            else
                _matchingAaSequences.Add(match, 1);
        }

        /// <summary>
        /// 在 peaklist 中查找 fi 离子，找到，返回其强度，否则返回 -1.
        /// </summary>
        /// <param name="sequence"></param>
        /// <param name="fi"></param>
        /// <param name="massTolerance"></param>
        /// <param name="peaklist"></param>
        /// <returns></returns>
        internal static double MatchPeak(AminoAcidSequence sequence, FragmentIon fi, double massTolerance, List<Peak> peaklist)
        {
            Peak peak = peaklist.Find(p => Math.Abs(fi.Mass - p.MassZ) <= massTolerance);
            if (peak != null && sequence != null)
                peak.addMatchingAASequence(sequence);
            return peak?.Intensity ?? -1.0;
        }

        internal static double MatchPeak(FragmentIon fi, double massTolerance, List<Peak> peaklist)
        {
            return MatchPeak(null, fi, massTolerance, peaklist);
        }

        // 2017-4-20
        public void ClearMatchData()
        {
            _matchingAaSequences = null;
        }

        public override string ToString()
        {
            return GetType().Name + "[mass=" + MassZ + ", intensity=" + Intensity + "]";
        }

        public int CompareTo(Peak p)
        {
            return MassZ.CompareTo(p.MassZ);
        }
    }
}