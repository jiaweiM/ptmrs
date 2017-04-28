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
    public class PTMResultClass : IDisposable
    {
        public struct Peptide
        {
            public int Id;
            public double Score;
            public double Probability;

            public Peptide(int id)
            {
                Id = id;
                Score = 0.0;
                Probability = 0.0;
            }

            public Peptide(int id, double score, double probability)
            {
                Id = id;
                Score = score;
                Probability = probability;
            }

            public override string ToString()
            {
                return $"{Id} {Score} {Probability}";
            }
        }

        public class IsoformGroup
        {
            public int SequenceID;

            public int SpectrumID;

            public readonly string SequenceString;

            public readonly string GroupReport;

            private Dictionary<AminoAcidModification, string> m_SiteProbabilityString;

            public bool Error { get; private set; }

            public string Message { get; private set; }

            public List<Peptide> Peptides { get; private set; }

            public List<int> PeptideIDs { get; private set; }

            public string BestSiteProbabilityString { get; private set; }

            public double Count
            {
                get
                {
                    double result = Error ? PeptideIDs.Count : Peptides.Count;
                    return result;
                }
            }

            public string Get_SiteProbabilityString(AminoAcidModification modification)
            {
                if (m_SiteProbabilityString == null)
                {
                    return null;
                }
                string result;
                if (m_SiteProbabilityString.TryGetValue(modification, out result))
                {
                    return result;
                }
                return null;
            }

            public IsoformGroup(string msg, List<int> peptideIDs)
            {
                Error = true;
                Message = msg;
                PeptideIDs = peptideIDs;
                Peptides = null;
            }

            public IsoformGroup(string sequenceString, Dictionary<AminoAcidModification, string> siteProbabilities,
                string bestSiteProbabilityString, string groupReport)
            {
                Error = false;
                Message = null;
                Peptides = new List<Peptide>();
                PeptideIDs = null;
                SequenceString = sequenceString;
                GroupReport = groupReport;
                m_SiteProbabilityString = siteProbabilities;
                BestSiteProbabilityString = bestSiteProbabilityString;
            }

            public void add(Peptide peptide)
            {
                Peptides.Add(peptide);
            }

            public override string ToString()
            {
                if (Error)
                {
                    return $"{SpectrumID} {Message}";
                }
                return $"{SpectrumID} {SequenceString} {GroupReport}";
            }
        }

        public delegate void ResultConsumer(
            BlockingCollection<IsoformGroup> result, CancellationTokenSource cl);

        private long _numberOfItems;

        internal BlockingCollection<IsoformGroup> IsoformGroupList;

        private bool _disposed;

        public long Count
        {
            get
            {
                long result;
                try
                {
                    result = Interlocked.Read(ref _numberOfItems);
                }
                catch (Exception)
                {
                    result = 0L;
                }
                return result;
            }
        }

        public PTMResultClass()
        {
            IsoformGroupList = new BlockingCollection<IsoformGroup>();
            _numberOfItems = 0L;
        }

        internal void CompleteAdding()
        {
            IsoformGroupList.CompleteAdding();
        }

        internal void AddRange(List<IsoformGroup> value, CancellationToken cl)
        {
            if (value != null)
            {
                foreach (IsoformGroup current in value)
                {
                    if (IsoformGroupList.TryAdd(current, -1, cl))
                    {
                        Interlocked.Increment(ref _numberOfItems);
                    }
                }
            }
        }

        internal void Add(IsoformGroup value, CancellationToken cl)
        {
            if (IsoformGroupList.TryAdd(value, -1, cl))
            {
                Interlocked.Increment(ref _numberOfItems);
            }
        }

        public static IDictionary<int, double> get_PeptideIdPrsScoreMap(List<IsoformGroup> groups)
        {
            IDictionary<int, double> result;
            int num = groups.Sum(delegate(IsoformGroup w)
            {
                if (!w.Error)
                {
                    return w.Peptides.Count;
                }
                return 0;
            });
            if (num == 0)
            {
                result = null;
            }
            else
            {
                Dictionary<int, double> dictionary = new Dictionary<int, double>(num);
                foreach (IsoformGroup current in from w in groups
                    where !w.Error
                    select w)
                {
                    foreach (Peptide current2 in current.Peptides)
                    {
                        dictionary.Add(current2.Id, current2.Score);
                    }
                }
                result = dictionary;
            }
            return result;
        }

        public static IDictionary<int, double> get_PeptideIdPrsProbabilityMap(List<IsoformGroup> groups)
        {
            IDictionary<int, double> result;
            int num = groups.Sum(delegate(IsoformGroup w)
            {
                if (!w.Error)
                {
                    return w.Peptides.Count;
                }
                return 0;
            });
            if (num == 0)
            {
                result = null;
            }
            else
            {
                Dictionary<int, double> dictionary = new Dictionary<int, double>(num);
                foreach (IsoformGroup current in from w in groups
                    where !w.Error
                    select w)
                {
                    foreach (Peptide current2 in current.Peptides)
                    {
                        dictionary.Add(current2.Id, current2.Probability);
                    }
                }
                result = dictionary;
            }
            return result;
        }

        public static IDictionary<int, string> get_PeptideIdPrsIsoformGroupReport(
            List<IsoformGroup> groups)
        {
            IDictionary<int, string> result;
            int num = groups.Sum(delegate(IsoformGroup w)
            {
                if (!w.Error)
                {
                    return w.Peptides.Count;
                }
                return 0;
            });
            if (num == 0)
            {
                result = null;
            }
            else
            {
                Dictionary<int, string> dictionary = new Dictionary<int, string>(num);
                foreach (IsoformGroup current in from w in groups
                    where !w.Error && !string.IsNullOrEmpty(w.GroupReport)
                    select w)
                {
                    foreach (Peptide current2 in current.Peptides)
                    {
                        dictionary.Add(current2.Id, current.GroupReport);
                    }
                }
                result = dictionary;
            }
            return result;
        }

        public static IDictionary<int, string> get_PeptideIdPrsSiteProbabilitiesMap(
            List<IsoformGroup> groups, AminoAcidModification aaModificiation)
        {
            IDictionary<int, string> result;
            int num = groups.Sum(delegate(IsoformGroup w)
            {
                if (!w.Error)
                {
                    return w.Peptides.Count;
                }
                return w.PeptideIDs.Count;
            });
            if (num == 0)
            {
                result = null;
            }
            else
            {
                Dictionary<int, string> dictionary = new Dictionary<int, string>(num);
                foreach (IsoformGroup current in groups)
                {
                    if (current.Error)
                    {
                        if (string.IsNullOrEmpty(current.Message))
                        {
                            continue;
                        }
                        using (List<int>.Enumerator enumerator2 = current.PeptideIDs.GetEnumerator())
                        {
                            while (enumerator2.MoveNext())
                            {
                                int current2 = enumerator2.Current;
                                dictionary.Add(current2, current.Message);
                            }
                            continue;
                        }
                    }
                    string value = current.Get_SiteProbabilityString(aaModificiation);
                    if (!string.IsNullOrEmpty(value))
                    {
                        foreach (int current3 in from w in current.Peptides
                            select w.Id)
                        {
                            dictionary.Add(current3, value);
                        }
                    }
                }
                result = dictionary;
            }
            return result;
        }

        public static IDictionary<int, string> get_PeptideIDBestPrsSiteProbabilityMap(
            List<IsoformGroup> groups)
        {
            IDictionary<int, string> result;
            int num = groups.Sum(delegate(IsoformGroup w)
            {
                if (!w.Error)
                {
                    return w.Peptides.Count;
                }
                return w.PeptideIDs.Count;
            });
            if (num == 0)
            {
                result = null;
            }
            else
            {
                Dictionary<int, string> dictionary = new Dictionary<int, string>(num);
                foreach (IsoformGroup current in groups)
                {
                    if (current.Error)
                    {
                        if (string.IsNullOrEmpty(current.Message))
                        {
                            continue;
                        }
                        using (List<int>.Enumerator enumerator2 = current.PeptideIDs.GetEnumerator())
                        {
                            while (enumerator2.MoveNext())
                            {
                                int current2 = enumerator2.Current;
                                dictionary.Add(current2, current.Message);
                            }
                            continue;
                        }
                    }
                    string bestSiteProbabilityString = current.BestSiteProbabilityString;
                    if (!string.IsNullOrEmpty(bestSiteProbabilityString))
                    {
                        foreach (int current3 in from w in current.Peptides
                            select w.Id)
                        {
                            dictionary.Add(current3, bestSiteProbabilityString);
                        }
                    }
                }
                result = dictionary;
            }
            return result;
        }

        public void Dispose()
        {
            Dispose(true);
            GC.SuppressFinalize(this);
        }

        protected virtual void Dispose(bool disposing)
        {
            if (!_disposed)
            {
                if (disposing && IsoformGroupList != null)
                {
                    IsoformGroupList.Dispose();
                    IsoformGroupList = null;
                }
                _disposed = true;
            }
        }
    }
}