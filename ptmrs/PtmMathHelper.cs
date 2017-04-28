#region License

// Created at 2016.07.13
// Author:  JiaweiMao
// Copyright (c) Dalian Institute of Chemical Physics
//            Chinese Academy of Sciences
// Contact: jiawei@dicp.ac.cn

#endregion

using System;
using System.Collections.Concurrent;
using System.Collections.Generic;
using System.Linq;
using System.Threading;

namespace ptmrs
{
    internal static class PtmMathHelper
    {
        private struct Snapshot
        {
            public int Elements;

            public int CurrBoxIndex;

            public int[] CurrBoxState;
        }

        private static readonly ConcurrentDictionary<long, double> Puffer = new ConcurrentDictionary<long, double>();

        private static long _numberOfItems;

        internal static int Count(this FitComposition value)
        {
            uint num1 = (uint)value;
            uint num = num1 - (num1 >> 1 & 1431655765U);
            num = (num & 858993459U) + (num >> 2 & 858993459u);
            return (int) ((num + (num >> 4) & 252645135u) * 16843009u >> 24);
        }

        /// <summary>
        /// Compare two double values.
        /// </summary>
        /// <returns>return 1 if a > b</returns>
        public static int Compare(double a, double b)
        {
            double num = a - b;
            if (num < 0.0)
                return -1;
            if (num <= 0.0)
                return 0;
            return 1;
        }

        // 4.20
        public static bool Equal(double a, double b, double accuracy)
        {
            return Math.Abs(a - b) < accuracy;
        }

        public static double GetPPMError(double trueValue, double measuredValue)
        {
            return (measuredValue - trueValue) / trueValue * 1000000.0;
        }

        public static double BinomialCoefficient(int n, int k)
        {
            if (k > n || n <= 0)
            {
                return 0.0;
            }
            if (k == 0)
            {
                return 1.0;
            }
            if (2 * k > n)
            {
                k = n - k;
            }
            long key = n * 2147483647L + k;
            double num;
            if (!Puffer.TryGetValue(key, out num))
            {
                num = 1.0;
                for (int i = 1; i <= k; i++)
                {
                    num *= ((double) n - k + i) / i;
                }
                if (Interlocked.Read(ref _numberOfItems) <= 10000L &&
                    Puffer.TryAdd(key, num))
                {
                    Interlocked.Increment(ref _numberOfItems);
                }
            }
            return num;
        }

        public static double CalculateBinomialScore(double probabilityRandamMatch, int numberTheoreticalFragments,
            int numberMatchedFragments)
        {
            if (numberMatchedFragments == 0 || probabilityRandamMatch == 0.0)
                return 0.0;
            double num = Math.Log(probabilityRandamMatch);
            double num2 = Math.Log(1.0 - probabilityRandamMatch);
            double probabilitySum = 0.0;
            for (int i = numberMatchedFragments; i <= numberTheoreticalFragments; i++)
            {
                probabilitySum += Math.Exp(Math.Log(BinomialCoefficient(numberTheoreticalFragments, i)) + i * num +
                             (numberTheoreticalFragments - i) * num2);
            }
            return Math.Max(0.0, -10.0 * Math.Log10(probabilitySum));
        }

        public static bool GetPossibleCombination(int n, int k, int[] maxN, out List<int[]> posPositions)
        {
            posPositions = new List<int[]>();
            if (n == 0)
            {
                return true;
            }
            if (k == 0)
            {
                return true;
            }
            if (maxN.Sum() == 0)
            {
                return true;
            }
            if (n > 15)
            {
                return false;
            }
            int num = maxN.Sum();
            int[] array = (from w in maxN
                select w + 1).ToArray();
            int num2 = 1;
            int[] array2 = array;
            for (int i = 0; i < array2.Length; i++)
            {
                int num3 = array2[i];
                num2 *= num3;
            }
            if (k > num)
            {
                return true;
            }
            if (k == num)
            {
                posPositions.Add(maxN);
                return true;
            }
            ConcurrentBag<int[]> concurrentBag = new ConcurrentBag<int[]>();
            for (int j = 0; j < num2; j++)
            {
                int[] item;
                if (CalculatePossibilities(j, array, k, out item))
                {
                    concurrentBag.Add(item);
                }
            }
            posPositions = concurrentBag.ToList();
            return true;
        }

        public static bool CalculatePossibilities(int id, int[] lengths, int numberOfModifications, out int[] indices)
        {
            int num = lengths.GetUpperBound(0) + 1;
            indices = new int[num];
            int num2 = 0;
            for (int i = num - 1; i >= 0; i--)
            {
                int num3 = 1;
                for (int j = 0; j < i; j++)
                {
                    num3 *= lengths[j];
                }
                int num4 = id % num3;
                indices[i] = (id - num4) / num3;
                id = num4;
                num2 += indices[i];
                if (num2 > numberOfModifications)
                {
                    indices = null;
                    return false;
                }
            }
            return num2 == numberOfModifications;
        }

        public static int[][] GetMultiCombinations(int elements, int boxes, int[] maxN, int maxNumberCombinations)
        {
            int num = maxN.Sum();
            if (elements > num)
            {
                return null;
            }
            if (elements == num)
            {
                return new[]
                {
                    maxN
                };
            }
            long num2 = 0L;
            List<int[][]> list = new List<int[][]>(10);
            int[][] array = new int[500][];
            Stack<Snapshot> stack = new Stack<Snapshot>();
            stack.Push(new Snapshot
            {
                CurrBoxIndex = 0,
                Elements = elements,
                CurrBoxState = new int[boxes]
            });
            while (stack.Count > 0)
            {
                Snapshot snapshot = stack.Pop();
                if (snapshot.CurrBoxIndex == boxes - 1)
                {
                    snapshot.CurrBoxState[snapshot.CurrBoxIndex] = snapshot.Elements;
                    if (num2 >= 500L)
                    {
                        list.Add(array);
                        array = new int[500][];
                        array[0] = snapshot.CurrBoxState;
                        num2 = 0L;
                    }
                    else
                    {
                        array[(int) (IntPtr) num2] = snapshot.CurrBoxState;
                    }
                    num2 += 1L;
                }
                else
                {
                    int num3 = 0;
                    for (int i = snapshot.CurrBoxIndex + 1; i < boxes; i++)
                    {
                        num3 += maxN[i];
                    }
                    int num4 = Math.Min(maxN[snapshot.CurrBoxIndex], snapshot.Elements);
                    for (int j = 0; j <= num4; j++)
                    {
                        if (snapshot.Elements - j <= num3)
                        {
                            Snapshot item = default(Snapshot);
                            if (snapshot.Elements - j > 0)
                            {
                                item.CurrBoxIndex = snapshot.CurrBoxIndex + 1;
                                item.Elements = snapshot.Elements - j;
                                if (j < num4)
                                {
                                    item.CurrBoxState = new int[boxes];
                                    Array.Copy(snapshot.CurrBoxState, item.CurrBoxState, boxes);
                                }
                                else
                                {
                                    item.CurrBoxState = snapshot.CurrBoxState;
                                }
                                item.CurrBoxState[snapshot.CurrBoxIndex] = j;
                            }
                            else
                            {
                                item.CurrBoxIndex = boxes - 1;
                                item.Elements = 0;
                                item.CurrBoxState = snapshot.CurrBoxState;
                                item.CurrBoxState[snapshot.CurrBoxIndex] = j;
                            }
                            stack.Push(item);
                        }
                    }
                }
            }
            if (num2 < 500L && num2 > 0L)
            {
                int[][] sourceArray = array;
                array = new int[num2][];
                Array.Copy(sourceArray, array, num2);
                list.Add(array);
            }
            num2 = 0L;
            for (int k = 0; k < list.Count; k++)
            {
                num2 += list[k].LongLength;
            }
            if (num2 > maxNumberCombinations)
            {
                return null;
            }
            int[][] array2 = new int[num2][];
            num2 = 0L;
            for (int l = 0; l < list.Count; l++)
            {
                for (int m = 0; m < list[l].Length; m++)
                {
                    array2[(int) ((IntPtr) num2)] = list[l][m];
                    num2 += 1L;
                }
            }
            return array2;
        }
    }
}