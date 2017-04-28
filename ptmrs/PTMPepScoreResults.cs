using System;
using System.Collections.Generic;
using System.Linq;

// 2017-04-19
namespace ptmrs
{
    /// <summary>
    /// PTM Peptide score for each window of a Peptide
    /// </summary>
    public class PTMPepScoreResults : IComparable<PTMPepScoreResults>
    {
        private readonly Dictionary<int, double> _scoreResults = new Dictionary<int, double>();

        /// <summary>
        /// The peak depth being compared.
        /// </summary>
        public int ComparisonPeakDepth { get; set; } = 1;

        /// <summary>
        /// 肽段序列
        /// </summary>
        public AminoAcidSequence Sequence { get; set; }

        public PTMPepScoreResults(AminoAcidSequence seq)
        {
            Sequence = seq;
        }

        /// <summary>
        /// 在指定 peak depth 对肽段按照打分降序排列
        /// </summary>
        public static void SortPTMPeptideScores(List<PTMPepScoreResults> pepScoreRes, int comparisonPeakDepth)
        {
            for (int i = 0; i < pepScoreRes.Count; i++)
                pepScoreRes[i].ComparisonPeakDepth = comparisonPeakDepth;
            pepScoreRes.Sort();
        }

        public override string ToString()
        {
            if (Sequence == null)
                return base.ToString();
            if (_scoreResults == null)
                return $"{Sequence} == NULL";
            return
                $"{Sequence} == {string.Join(",", from w in _scoreResults select w.ToString())}";
        }

        /// <summary>
        /// Add a peak depth and its peptide score value.
        /// </summary>
        public void AddScoreResult(int peakDepth, double score)
        {
            _scoreResults.Add(peakDepth, score);
        }

        /// <summary>
        /// 获得指定 peak depth 下的打分值。
        /// </summary>
        public double GetScoreResult(int peakDepth)
        {
            double result;
            if (_scoreResults.TryGetValue(peakDepth, out result))
                return result;
            return -1.0;
        }

        /// <summary>
        /// 获得最高的打分值及对应的 peak depth. 4.20
        /// </summary>
        /// <param name="minPeakDepth">输出的 peak depth 的默认值</param>
        /// <param name="outPeakDepth">输出的最高打分值的 peak depth</param>
        /// <returns>最高打分值</returns>
        public double GetMaxScoreResult(int minPeakDepth, out int outPeakDepth)
        {
            if (_scoreResults == null || _scoreResults.Count == 0)
            {
                outPeakDepth = -1;
                return -1.0;
            }
            double score = 0.0;
            outPeakDepth = minPeakDepth;
            foreach (KeyValuePair<int, double> current in _scoreResults)
            {
                if (current.Value > score)
                {
                    score = current.Value;
                    outPeakDepth = current.Key;
                }
            }
            return score;
        }

        public int CompareTo(PTMPepScoreResults scoreRes)
        {
            if (_scoreResults == null || _scoreResults.Count == 0 || scoreRes._scoreResults == null ||
                scoreRes._scoreResults.Count == 0)
                return 0;
            return CompareTo(scoreRes, ComparisonPeakDepth);
        }

        /// <summary>
        /// Compare with another PTMPepScoreResults with given peak depth.
        /// </summary>
        /// <returns>如果当前值比比较值大，返回1，相等返回0，否则返回-1</returns>
        public int CompareTo(PTMPepScoreResults score, int peakDepth)
        {
            double b;
            double a;
            if (_scoreResults.TryGetValue(peakDepth, out b) && score._scoreResults.TryGetValue(peakDepth, out a))
                return PtmMathHelper.Compare(a, b);
            return 0;
        }
    }
}