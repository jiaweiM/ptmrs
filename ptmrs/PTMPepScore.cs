#region License

// Created at 2016.07.13
// Author:  JiaweiMao
// Copyright (c) Dalian Institute of Chemical Physics
//            Chinese Academy of Sciences
// Contact: jiawei@dicp.ac.cn

#endregion

using System;
using System.Collections.Generic;
using System.Linq;

namespace ptmrs
{
    /// <summary>
    /// 
    /// </summary>
    public class PTMPepScore : IComparable<PTMPepScore>
    {
        public AminoAcidSequence Sequence { get; private set; }
        public double Score { get; private set; }
        public double ReciprocalScore { get; private set; }

        public PTMPepScore(AminoAcidSequence seq, double score)
        {
            Sequence = seq;
            Score = score;
            ReciprocalScore = Math.Pow(10.0, Score / 10.0);
        }

        public override string ToString()
        {
            if (Sequence == null)
            {
                return base.ToString();
            }
            return $"{Sequence}, Score: {Score}";
        }

        public int CompareTo(PTMPepScore o)
        {
            return PtmMathHelper.Compare(o.Score, Score);
        }

        /// <summary>
        /// Return the sum of reciprocal score.
        /// </summary>
        /// <param name="ptmPepScores">List of PTMPepScore.</param>
        /// <returns>sum of reciprocal score.</returns>
        public static double GetReciprocalPSum(List<PTMPepScore> ptmPepScores)
        {
            return ptmPepScores?.Sum(w => w.ReciprocalScore) ?? 0.0;
        }
    }
}