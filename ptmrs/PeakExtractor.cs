using System;
using System.Collections.Generic;
using System.Linq;

namespace ptmrs
{
    public class PeakExtractor
    {
        // 4.19
        public class MZWindow
        {
            public enum PeakAddedResult
            {
                Added,
                DiscardedMzLowerLimit,
                DiscardedMzUpperLimit
            }

            private bool _isSorted;
            private readonly int _lastPeakExtractIndex;

            public Peak this[int n] => GetNMostIntensePeak(n);

            public double LowerLimit { get; }
            public double UpperLimit { get; }

            public List<Peak> Peaks { get; } = new List<Peak>();

            public MZWindow(double lowerLimit, double upperLimit)
            {
                _isSorted = false;
                _lastPeakExtractIndex = -1;
                LowerLimit = lowerLimit;
                UpperLimit = upperLimit;
            }

            public bool HasNextMostIntensePeak()
            {
                return !IsEmpty() && _lastPeakExtractIndex + 1 < Peaks.Count;
            }

            public PeakAddedResult AddPeak(Peak p)
            {
                if (p.MassZ < LowerLimit)
                    return PeakAddedResult.DiscardedMzLowerLimit;
                if (p.MassZ >= UpperLimit)
                    return PeakAddedResult.DiscardedMzUpperLimit;
                Peaks.Add(p);
                _isSorted = false;
                return PeakAddedResult.Added;
            }

            private void SortPeaksDescendingByIntensity()
            {
                Peaks.SortDescendingByIntensity();
                _isSorted = true;
            }

            /// <summary>
            /// 获得强度前 peakDepth 的谱峰
            /// </summary>
            public List<Peak> GetMostIntensePeaks(int peakDepth)
            {
                if (Peaks == null || peakDepth > Peaks.Count)
                    return Peaks;
                if (!_isSorted)
                    SortPeaksDescendingByIntensity();
                return Peaks.GetRange(0, peakDepth);
            }

            /// <summary>
            /// 获得强度排第 n 的谱峰 
            /// </summary>
            private Peak GetNMostIntensePeak(int n)
            {
                if (Peaks == null || n > Peaks.Count)
                    return null;
                if (!_isSorted)
                    SortPeaksDescendingByIntensity();
                return Peaks[n - 1];
            }

            /// <summary>
            /// 获得最小强度及其 index.
            /// </summary>
            public double GetMinIntensity(out int outMinIntensityIndex)
            {
                double num = int.MaxValue;
                outMinIntensityIndex = -1;
                for (int i = 0; i < Peaks.Count; i++)
                {
                    Peak peak = Peaks[i];
                    if (peak.Intensity < num)
                    {
                        num = peak.Intensity;
                        outMinIntensityIndex = i;
                    }
                }
                return num;
            }

            public bool IsEmpty()
            {
                return Peaks == null || Peaks.Count == 0;
            }
        }

        protected List<MZWindow> MZWindows;
        public double WindowWidth { get; }
        public double StartMassZ { get; }
        public double EndMassZ { get; }
        public int MzWindowCount => MZWindows.Count;

        public PeakExtractor(PeptideSpectrumMatch spec, int minPeakDepth, int maxPeakDepth)
        {
            WindowWidth = 100.0;
            double startMassZ = 50.0;
            double startMassZ2;
            double endMassZ;
            List<Peak> peaks = spec.Peaks;
            MZWindows = DistributePeaksToMassZWindows(peaks, startMassZ, WindowWidth,
                out startMassZ2, out endMassZ);
            StartMassZ = startMassZ2;
            EndMassZ = endMassZ;
        }

        private List<Peak> ExtractMostIntensePeaks(int window, int peakDepth)
        {
            if (MZWindows == null || window < 0 || window >= MZWindows.Count)
            {
                return new List<Peak>(1);
            }
            return MZWindows[window].GetMostIntensePeaks(peakDepth);
        }

        public Peak GetNMostIntensePeaks(int window, int nMostIntense)
        {
            return MZWindows[window][nMostIntense];
        }

        public List<Peak> GetMostIntensePeaks(int peakDepth)
        {
            return
            (from p in
                MZWindows.SelectMany(wnd => wnd.GetMostIntensePeaks(peakDepth))
                orderby p.Intensity descending
                select p).ToList();
        }

        public List<Peak> GetMostIntensePeaks(List<int> peakDepths)
        {
            List<Peak> list = new List<Peak>(peakDepths.Sum());
            for (int i = 0; i < peakDepths.Count; i++)
            {
                list.AddRange(MZWindows[i].GetMostIntensePeaks(peakDepths[i]));
            }
            return list;
        }

        public Tuple<double, double> GetMzWindowBounds(int window)
        {
            if (MZWindows == null || window < 0 || window >= MZWindows.Count)
            {
                return null;
            }
            return Tuple.Create(MZWindows[window].LowerLimit, MZWindows[window].UpperLimit);
        }

        public List<MZWindow> DistributePeaksToMassZWindows<TP>(List<TP> rawPeaks, double startMassZ,
            double windowWidth, out double outStartMassZ, out double outEndMassZ) where TP : Peak
        {
            TP p = rawPeaks[0];
            TP p2 = rawPeaks[rawPeaks.Count - 1];
            int num = 0;
            if (p.MassZ > startMassZ)
            {
                num = (int) ((p.MassZ - startMassZ) / windowWidth);
            }
            List<MZWindow> list =
                new List<MZWindow>((int) ((p2.MassZ - startMassZ) / windowWidth) + 1 - num);
            MZWindow mZWindow = new MZWindow(startMassZ + windowWidth * num,
                startMassZ + windowWidth * (num + 1));
            list.Add(mZWindow);
            foreach (TP current in rawPeaks)
            {
                MZWindow.PeakAddedResult peakAddedResult = mZWindow.AddPeak(current);
                if (peakAddedResult == MZWindow.PeakAddedResult.DiscardedMzUpperLimit)
                {
                    do
                    {
                        mZWindow = new MZWindow(mZWindow.UpperLimit, mZWindow.UpperLimit + windowWidth);
                        list.Add(mZWindow);
                        peakAddedResult = mZWindow.AddPeak(current);
                    } while (peakAddedResult == MZWindow.PeakAddedResult.DiscardedMzUpperLimit);
                }
            }
            outStartMassZ = list[0].LowerLimit;
            outEndMassZ = list[list.Count - 1].UpperLimit;
            return list;
        }
    }
}