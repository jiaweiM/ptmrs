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
    internal static class PtmMathExtensioEnumerableExtensions
    {
        private class AllPossibleCombination
        {
            private int m_n;

            private int m_k;

            private int[] m_indices;

            private List<int[]> m_combinations;

            public AllPossibleCombination(int n, int k)
            {
                if (n <= 0)
                {
                    throw new ArgumentException("n_ must be in N+");
                }
                if (k <= 0)
                {
                    throw new ArgumentException("k_ must be in N+");
                }
                if (k > n)
                {
                    throw new ArgumentException("k_ can be at most n_");
                }
                this.m_n = n;
                this.m_k = k;
                this.m_indices = new int[this.m_k];
                this.m_indices[0] = 1;
            }

            public List<int[]> GetCombinations()
            {
                if (this.m_combinations == null)
                {
                    this.m_combinations = new List<int[]>();
                    this.Iterate(0);
                }
                return this.m_combinations;
            }

            private void Iterate(int ii)
            {
                if (ii > 0)
                {
                    this.m_indices[ii] = this.m_indices[ii - 1] + 1;
                }
                while (this.m_indices[ii] <= this.m_n - this.m_k + ii + 1)
                {
                    if (ii < this.m_k - 1)
                    {
                        this.Iterate(ii + 1);
                    }
                    else
                    {
                        int[] array = new int[this.m_k];
                        this.m_indices.CopyTo(array, 0);
                        this.m_combinations.Add(array);
                    }
                    this.m_indices[ii]++;
                }
            }
        }

        public class Combinations
        {
            public List<List<Tuple<AminoAcidModification, List<ModificationPosition>>>> possibleCombination;

            public Combinations(Dictionary<AminoAcidModification, List<List<ModificationPosition>>> values)
            {
                List<List<Tuple<AminoAcidModification, List<ModificationPosition>>>> list = new List<List<Tuple<AminoAcidModification, List<ModificationPosition>>>>(values.Count);
                int num = 0;
                foreach (KeyValuePair<AminoAcidModification, List<List<ModificationPosition>>> current in values)
                {
                    List<Tuple<AminoAcidModification, List<ModificationPosition>>> list2 = new List<Tuple<AminoAcidModification, List<ModificationPosition>>>(current.Value.Count);
                    for (int i = 0; i < current.Value.Count; i++)
                    {
                        Tuple<AminoAcidModification, List<ModificationPosition>> item = new Tuple<AminoAcidModification, List<ModificationPosition>>(current.Key, current.Value[i]);
                        list2.Add(item);
                    }
                    list.Add(list2);
                    num++;
                }
                this.calculate(list);
            }

            private void calculate(List<List<Tuple<AminoAcidModification, List<ModificationPosition>>>> values)
            {
                try
                {
                    List<List<Tuple<AminoAcidModification, List<ModificationPosition>>>> list = PtmMathExtensioEnumerableExtensions.Combinations.AllCombinationsOf<Tuple<AminoAcidModification, List<ModificationPosition>>>(values.ToArray());
                    this.possibleCombination = new List<List<Tuple<AminoAcidModification, List<ModificationPosition>>>>(list.Count);
                    for (int i = 0; i < list.Count; i++)
                    {
                        bool flag = false;
                        for (int j = 0; j < list[i].Count; j++)
                        {
                            for (int k = j + 1; k < list[i].Count; k++)
                            {
                                if (!this.comparer(list[i][j].Item2, list[i][k].Item2))
                                {
                                    flag = true;
                                    break;
                                }
                            }
                        }
                        if (!flag)
                        {
                            this.possibleCombination.Add(list[i]);
                        }
                    }
                    this.possibleCombination.TrimExcess();
                }
                catch (Exception)
                {
                    throw;
                }
            }

            private bool comparer(List<ModificationPosition> value1, List<ModificationPosition> value2)
            {
                for (int i = 0; i < value1.Count; i++)
                {
                    for (int j = 0; j < value2.Count; j++)
                    {
                        if (value1[i].ZeroBasedPosition == value2[j].ZeroBasedPosition)
                        {
                            return false;
                        }
                    }
                }
                return true;
            }

            private static List<List<T>> AllCombinationsOf<T>(params List<T>[] sets)
            {
                List<List<T>> list = new List<List<T>>();
                foreach (T current in sets[0])
                {
                    list.Add(new List<T>
                    {
                        current
                    });
                }
                foreach (List<T> current2 in sets.Skip(1))
                {
                    list = PtmMathExtensioEnumerableExtensions.Combinations.AddExtraSet<T>(list, current2);
                }
                return list;
            }

            private static List<List<T>> AddExtraSet<T>(List<List<T>> combinations, List<T> set)
            {
                IEnumerable<List<T>> source = from value in set
                                              from combination in combinations
                                              select new List<T>(combination)
                {
                    value
                };
                return source.ToList<List<T>>();
            }
        }

        public static List<List<E>> GetCombinationsWithoutReplication<E>(this List<E> list, int elementsToChoose)
        {
            List<int[]> combinations = new PtmMathExtensioEnumerableExtensions.AllPossibleCombination(list.Count, elementsToChoose).GetCombinations();
            List<List<E>> list2 = new List<List<E>>();
            foreach (int[] current in combinations)
            {
                List<E> list3 = new List<E>(current.Length);
                int[] array = current;
                for (int i = 0; i < array.Length; i++)
                {
                    int num = array[i];
                    list3.Add(list[num - 1]);
                }
                list2.Add(list3);
            }
            return list2;
        }
    }
}