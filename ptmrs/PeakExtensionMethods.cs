using System.Collections.Generic;

// 2017-04-19
namespace ptmrs
{
    public static class PeakExtensionMethods
    {
        /// <summary>
        /// 谱峰按强度降序排列。
        /// </summary>
        /// <param name="peaks">待排序谱峰</param>
        public static void SortDescendingByIntensity(this List<Peak> peaks)
        {
            peaks.Sort((peak1, peak2) => peak2.Intensity.CompareTo(peak1.Intensity));
        }

        /// <summary>
        /// 查找特定的谱峰，即 m/z 范围内强度最高的谱峰
        /// </summary>
        /// <param name="mass">待查找谱峰 m/z</param>
        /// <param name="deltaMass">Tolerance窗口宽度</param>
        /// <returns></returns>
        public static TE Find<TE>(this IEnumerable<TE> peaks, double mass, double deltaMass) where TE : Peak
        {
            bool flag = false;
            TE e2 = default(TE);
            foreach (TE current in peaks)
            {
                if (mass - deltaMass / 2.0 <= current.MassZ && current.MassZ <= mass + deltaMass / 2.0 &&
                    (e2 == null || current.Intensity > e2.Intensity))
                {
                    flag = true;
                    e2 = current;
                }
            }
            if (!flag)
            {
                return default(TE);
            }
            return e2;
        }
    }
}