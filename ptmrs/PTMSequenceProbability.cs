#region License

// Created at 2016.07.13
// Author:  JiaweiMao
// Copyright (c) Dalian Institute of Chemical Physics
//            Chinese Academy of Sciences
// Contact: jiawei@dicp.ac.cn

#endregion

using System;

namespace ptmrs
{
    public class PTMSequenceProbability : IComparable, IComparable<PTMSequenceProbability>
    {
        public AminoAcidSequence Sequence { get; }
        public double Score { get; }
        public double Probability { get; }

        public PTMSequenceProbability(AminoAcidSequence seq, double pepScore, double probability)
        {
            Sequence = seq;
            Score = pepScore;
            Probability = probability;
        }

        public override string ToString()
        {
            return $"{Sequence}, Score: {Score}, Probability: {Probability}%";
        }

        public int CompareTo(object obj)
        {
            if (obj == null)
            {
                return -1;
            }
            PTMSequenceProbability pTmSequenceProbability = obj as PTMSequenceProbability;
            if (pTmSequenceProbability == null)
            {
                throw new ArgumentException("Argument must be PTMSequenceProbability");
            }
            return CompareTo(pTmSequenceProbability);
        }

        public int CompareTo(PTMSequenceProbability other)
        {
            if (other == null)
            {
                return -1;
            }
            return other.Score.CompareTo(Score);
        }
    }
}