using System;
using System.Collections.Concurrent;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text;
using System.Threading;
using System.Threading.Tasks;
using System.Xml.Linq;
using ptmrs.Data;
using ptmrs.Properties;

namespace ptmrs
{
    public class AnyPTMProcessor : IDisposable
    {
        private class PeptideSpectrumsSource : IDisposable
        {
            private const int m_boundarySizeOfPackage64bit = 3000;

            private const int m_boundarySizeOfPackage32bit = 1500;

            private bool disposed;

            private IDataConection m_DataConnection;

            private BlockingCollection<progressMessage> progressMessageQueue;

            private int m_searchSpectraUnder;

            private readonly int m_desiredSizeOfPackages;

            private long m_numberOfSpectraInPackage;

            private ConcurrentQueue<List<SpectraPackageItem>> m_PeptideSpectrums;

            private AutoResetEvent m_Barrier;

            private ManualResetEvent m_Manualreset_isCollectingSpectra;

            private ManualResetEvent m_Manualreset_spectraLeft;

            public int MaxTaskCount { get; private set; }

            private bool IsSpectraLeft => m_Manualreset_spectraLeft.WaitOne(0);

            private bool IsCollectingSpectra => m_Manualreset_isCollectingSpectra.WaitOne(0);

            public PeptideSpectrumsSource(int totoalNumberOfSpectra, IDataConection dataConnection,
                BlockingCollection<progressMessage> progressMessageQueue, int maxTaskCount)
            {
                m_desiredSizeOfPackages = 0;
                m_DataConnection = dataConnection;
                this.progressMessageQueue = progressMessageQueue;
                int num;
                CalculateOptimalSizeOfPackage(totoalNumberOfSpectra, out num, out m_desiredSizeOfPackages,
                    out m_searchSpectraUnder);
                if (maxTaskCount > 0)
                    MaxTaskCount = Math.Min(maxTaskCount, num);
                else
                    MaxTaskCount = num;
                progressMessageQueue.Add(
                    new progressMessage(string.Format("Workload level: #spectra: {0}. ", totoalNumberOfSpectra)));
                progressMessageQueue.Add(
                    new progressMessage(string.Format("Workload level: #spectra per package: {0}. ",
                        m_desiredSizeOfPackages)));
                progressMessageQueue.Add(
                    new progressMessage(string.Format("Workload level: #parallel tasks: {0}. ", MaxTaskCount)));
                m_Barrier = new AutoResetEvent(true);
                m_Manualreset_isCollectingSpectra = new ManualResetEvent(false);
                m_Manualreset_spectraLeft = new ManualResetEvent(true);
                m_PeptideSpectrums = new ConcurrentQueue<List<SpectraPackageItem>>();
            }

            private void GetNewPackage()
            {
                try
                {
                    m_Barrier.WaitOne();
                    long num = Interlocked.Read(ref m_numberOfSpectraInPackage);
                    if (IsSpectraLeft && num <= m_searchSpectraUnder)
                    {
                        try
                        {
                            m_Manualreset_isCollectingSpectra.Set();
                            int num2;
                            int num3;
                            List<List<SpectraPackageItem>> newDataPackage =
                                m_DataConnection.GetNewDataPackage(m_desiredSizeOfPackages, out num2,
                                    out num3);
                            if (newDataPackage == null)
                            {
                                progressMessageQueue.TryAdd(new progressMessage("Finished collecting spectra"));
                                m_Manualreset_spectraLeft.Reset();
                            }
                            else
                            {
                                Parallel.ForEach(newDataPackage,
                                    delegate(List<SpectraPackageItem> item) { m_PeptideSpectrums.Enqueue(item); });
                                Interlocked.Add(ref m_numberOfSpectraInPackage,
                                    m_PeptideSpectrums.Count);
                                m_Manualreset_spectraLeft.Set();
                            }
                        }
                        finally
                        {
                            m_Manualreset_isCollectingSpectra.Reset();
                        }
                    }
                }
                catch (Exception)
                {
                    m_Manualreset_spectraLeft.Reset();
                    throw;
                }
                finally
                {
                    m_Barrier.Set();
                }
            }

            private static void CalculateOptimalSizeOfPackage(int totalNumberOfSpectra, out int maxTaskCount,
                out int desiredSizeOfPackages, out int searchSpectraUnder)
            {
                maxTaskCount = Environment.ProcessorCount;
                int num;
                if (!Environment.Is64BitProcess)
                {
                    num = 1500;
                }
                else
                {
                    num = 3000;
                }
                desiredSizeOfPackages = Math.Max(Math.Min(num * maxTaskCount, totalNumberOfSpectra), 40000);
                searchSpectraUnder = (int) Math.Max(desiredSizeOfPackages * 0.3, 200.0);
            }

            public List<SpectraPackageItem> GetNext()
            {
                List<SpectraPackageItem> result;
                List<SpectraPackageItem> next;
                if (m_PeptideSpectrums.TryDequeue(out next))
                {
                    Interlocked.Decrement(ref m_numberOfSpectraInPackage);
                    progressMessageQueue.TryAdd(
                        new progressMessage(next.Sum(w => w.SpectraCount),
                            next.Sum(w => w.PeptideCount)));
                    if (!IsCollectingSpectra && Interlocked.Read(ref m_numberOfSpectraInPackage) <=
                        m_searchSpectraUnder)
                    {
                        GetNewPackage();
                    }
                }
                if (next == null)
                {
                    GetNewPackage();
                }
                if (IsSpectraLeft && next == null)
                {
                    next = GetNext();
                }
                result = next;
                return result;
            }

            public void Dispose()
            {
                Dispose(true);
                GC.SuppressFinalize(this);
            }

            protected virtual void Dispose(bool disposing)
            {
                if (!disposed)
                {
                    if (disposing)
                    {
                        if (m_Manualreset_isCollectingSpectra != null)
                        {
                            m_Manualreset_isCollectingSpectra.Dispose();
                            m_Manualreset_isCollectingSpectra = null;
                        }
                        if (m_Manualreset_spectraLeft != null)
                        {
                            m_Manualreset_spectraLeft.Dispose();
                            m_Manualreset_spectraLeft = null;
                        }
                        if (m_Barrier != null)
                        {
                            m_Barrier.Dispose();
                            m_Barrier = null;
                        }
                        if (m_PeptideSpectrums != null)
                        {
                            m_PeptideSpectrums = null;
                        }
                    }
                    disposed = true;
                }
            }
        }

        public PTMResultClass PtmResult;

        private bool disposed;

        private int m_MaximalPeakDepth;

        private int m_RandomSeed;

        private bool m_UseDiagnosticIons;

        private double m_MassTolerance;

        private int _maxPtmCount;

        private int _maxIsoformCount;

        private ThreeStateEnum m_ScoreNLPeaksToo;

        private CancellationTokenSource m_SearchCancel;

        private PeptideSpectrumsSource _psmSource;

        private readonly BlockingCollection<progressMessage> progressMessageQueue;

        private readonly ModificationManagement _modificationManagment;

        private Dictionary<SpectrumType, FitComposition> _nLsFitComposition;

        private Dictionary<SpectrumType, FitComposition> _fitCompositions;

        public bool ErrorHappened { get; private set; }

        public bool UseMassAccuracyCorrection { get; set; }

        public ModificationManagement EquivalentModificationMap => _modificationManagment;

        internal FitComposition getNLsFragmentIonComposition(SpectrumType activationType)
        {
            if (m_ScoreNLPeaksToo == ThreeStateEnum.False)
            {
                return FitComposition.No;
            }
            FitComposition result;
            if (_nLsFitComposition != null && _nLsFitComposition.TryGetValue(activationType, out result))
                return result;
            switch (activationType)
            {
                case SpectrumType.CID_CAD:
                    if (m_ScoreNLPeaksToo == ThreeStateEnum.True)
                        return FitComposition.B | FitComposition.Y;
                    return FitComposition.No;
                case SpectrumType.ECD_ETD:
                    return FitComposition.No;
                case SpectrumType.HCD:
                case SpectrumType.EThcD:
                    return FitComposition.B | FitComposition.Y;
            }
            throw new ArgumentException("Spectrum type unknown");
        }

        internal FitComposition GetFragmentIonTypeComposition(SpectrumType activationType)
        {
            FitComposition result;
            if (_fitCompositions != null && _fitCompositions.TryGetValue(activationType, out result))
                return result;
            switch (activationType)
            {
                case SpectrumType.None:
                    throw new ArgumentException("Spectrum type no set");
                case SpectrumType.CID_CAD:
                case SpectrumType.HCD:
                    return FitComposition.B | FitComposition.Y;
                case SpectrumType.ECD_ETD:
                    return FitComposition.C | FitComposition.ZRadical | FitComposition.ZPrime;
                case SpectrumType.EThcD:
                    return FitComposition.B | FitComposition.Y | FitComposition.C | FitComposition.ZRadical |
                           FitComposition.ZPrime;
                default:
                    throw new ArgumentException("Spectrum type unknown");
            }
        }

        private void ReadXMLFragmentIonTypeComposition(XDocument doc)
        {
            XElement xElement = doc.Element("FragmentIonCompositionPreference");
            if (xElement == null)
                return;

            _fitCompositions = new Dictionary<SpectrumType, FitComposition>();
            _nLsFitComposition = new Dictionary<SpectrumType, FitComposition>();
            foreach (XElement current in from w in xElement.Elements()
                where string.Compare(w.Name.LocalName, "FragmentIonComposition") == 0
                select w)
            {
                string value;
                if ((value = current.Attribute("ActivationType").Value) != null)
                {
                    if (< PrivateImplementationDetails >{
                        C2D0B23E - 2264 - 4CE1 - B95A - E19837CB9B1B
                    }.$$method0x60001dc - 1 == null)
                    {
                        < PrivateImplementationDetails >{
                            C2D0B23E - 2264 - 4CE1 - B95A - E19837CB9B1B
                        }.$$method0x60001dc - 1 = new Dictionary<string, int>(8)
                        {
                            {
                                "CID_CAD",
                                0
                            },
                            {
                                "CAD",
                                1
                            },
                            {
                                "CID",
                                2
                            },
                            {
                                "ETD",
                                3
                            },
                            {
                                "ECD",
                                4
                            },
                            {
                                "ECD_ETD",
                                5
                            },
                            {
                                "EThcD",
                                6
                            },
                            {
                                "HCD",
                                7
                            }
                        };
                    }
                    int num;
                    if (< PrivateImplementationDetails >{
                        C2D0B23E - 2264 - 4CE1 - B95A - E19837CB9B1B
                    }.$$method0x60001dc - 1.TryGetValue(value, out num))
                    {
                        SpectrumType spectrumType;
                        switch (num)
                        {
                            case 0:
                            case 1:
                            case 2:
                                spectrumType = SpectrumType.CID_CAD;
                                break;
                            case 3:
                            case 4:
                            case 5:
                                spectrumType = SpectrumType.ECD_ETD;
                                break;
                            case 6:
                                spectrumType = SpectrumType.EThcD;
                                break;
                            case 7:
                                spectrumType = SpectrumType.HCD;
                                break;
                            default:
                                continue;
                        }
                        FitComposition fITComposition = FitComposition.No;
                        string[] array = current.Attribute("FragmentIonComposition")
                            .Value.Split(',');
                        for (int i = 0; i < array.Length; i++)
                        {
                            string text = array[i];
                            string a;
                            if ((a = text) != null)
                            {
                                if (!(a == "b"))
                                {
                                    if (!(a == "y"))
                                    {
                                        if (!(a == "c"))
                                        {
                                            if (!(a == "zPrime"))
                                            {
                                                if (a == "zRadical")
                                                {
                                                    fITComposition |= FitComposition.ZRadical;
                                                }
                                            }
                                            else
                                            {
                                                fITComposition |= FitComposition.ZPrime;
                                            }
                                        }
                                        else
                                        {
                                            fITComposition |= FitComposition.C;
                                        }
                                    }
                                    else
                                    {
                                        fITComposition |= FitComposition.Y;
                                    }
                                }
                                else
                                {
                                    fITComposition |= FitComposition.B;
                                }
                            }
                        }
                        FitComposition fITComposition2 = FitComposition.No;
                        string[] array2 = current.Attribute("NeutralLossFragmentIonComposition")
                            .Value.Split(',');
                        for (int j = 0; j < array2.Length; j++)
                        {
                            string text2 = array2[j];
                            string a2;
                            if ((a2 = text2) != null)
                            {
                                if (!(a2 == "b"))
                                {
                                    if (!(a2 == "y"))
                                    {
                                        if (!(a2 == "c"))
                                        {
                                            if (!(a2 == "zPrime"))
                                            {
                                                if (a2 == "zRadical")
                                                {
                                                    fITComposition2 |= FitComposition.ZRadical;
                                                }
                                            }
                                            else
                                            {
                                                fITComposition2 |= FitComposition.ZPrime;
                                            }
                                        }
                                        else
                                        {
                                            fITComposition2 |= FitComposition.C;
                                        }
                                    }
                                    else
                                    {
                                        fITComposition2 |= FitComposition.Y;
                                    }
                                }
                                else
                                {
                                    fITComposition2 |= FitComposition.B;
                                }
                            }
                        }
                        FitComposition fragmentIonTypeComposition =
                            GetFragmentIonTypeComposition(spectrumType);
                        if (fragmentIonTypeComposition != fITComposition)
                        {
                            _fitCompositions[spectrumType] = fITComposition;
                        }
                        FitComposition nLsFragmentIonComposition =
                            getNLsFragmentIonComposition(spectrumType);
                        if (nLsFragmentIonComposition != fITComposition2)
                        {
                            _nLsFitComposition[spectrumType] = fITComposition2;
                        }
                    }
                }
            }
        }

        public static string DefaulConfigurationXML()
        {
            return Resources.IMP_ptmRSConf;
        }

        private void initialize(XDocument doc, IDataConection dataConnection, int totalNumberOfSpectra,
            int maxTaskCount, List<AminoAcidModification> scoredAAMods)
        {
            _modificationManagment.ReadXML_EquivalentModifications(doc, scoredAAMods);
            List<string> list;
            _modificationManagment.CompleteEquivalentModificationAdding(scoredAAMods, m_MassTolerance,
                out list);
            if (list != null)
            {
                foreach (string current in list)
                {
                    progressMessageQueue.TryAdd(new progressMessage(current));
                }
            }
            ReadXMLFragmentIonTypeComposition(doc);
            progressMessageQueue.TryAdd(new progressMessage(string.Format("FITs for {2}: {0}; FITs with NLs: {1}",
                GetFragmentIonTypeComposition(SpectrumType.CID_CAD),
                getNLsFragmentIonComposition(SpectrumType.CID_CAD), SpectrumType.CID_CAD)));
            progressMessageQueue.TryAdd(new progressMessage(string.Format("FITs for {2}: {0}; FITs with NLs: {1}",
                GetFragmentIonTypeComposition(SpectrumType.HCD),
                getNLsFragmentIonComposition(SpectrumType.HCD), SpectrumType.HCD)));
            progressMessageQueue.TryAdd(new progressMessage(string.Format("FITs for {2}: {0}; FITs with NLs: {1}",
                GetFragmentIonTypeComposition(SpectrumType.ECD_ETD),
                getNLsFragmentIonComposition(SpectrumType.ECD_ETD), SpectrumType.ECD_ETD)));
            progressMessageQueue.TryAdd(new progressMessage(string.Format("FITs for {2}: {0}; FITs with NLs: {1}",
                GetFragmentIonTypeComposition(SpectrumType.EThcD),
                getNLsFragmentIonComposition(SpectrumType.EThcD), SpectrumType.EThcD)));
            _psmSource = new PeptideSpectrumsSource(totalNumberOfSpectra, dataConnection,
                progressMessageQueue, maxTaskCount);
        }

        public AnyPTMProcessor(IDataConection dataConnection, CancellationTokenSource m_searchCancel,
            string configFileName, int RandomSeed, bool UseDiagnosticIons, int maxPeakDepth, int maxIsoformCount,
            int maxPTMCount, ThreeStateEnum scoreNLPeaksToo, double masstolerance,
            List<AminoAcidModification> scoredAAMods, int totalNumberOfSpectra, int maxTaskCount)
        {
            progressMessageQueue = dataConnection.GetProgressMessageQueue();
            if (progressMessageQueue == null)
                throw new ArgumentNullException(
                    "ProgressmessageQueue is null! Please provide a Blocking Collection through the function GetProgressMessageQueue.");
            m_MaximalPeakDepth = maxPeakDepth;
            m_RandomSeed = RandomSeed;
            m_UseDiagnosticIons = UseDiagnosticIons;
            _maxIsoformCount = maxIsoformCount;
            _maxPtmCount = maxPTMCount;
            m_ScoreNLPeaksToo = scoreNLPeaksToo;
            m_MassTolerance = masstolerance;
            ErrorHappened = false;
            m_SearchCancel = m_searchCancel;
            PtmResult = new PTMResultClass();
            _modificationManagment = new ModificationManagement();
            FileInfo fileInfo = new FileInfo(configFileName);
            if (!fileInfo.Exists)
                throw new FileNotFoundException("Configuration file not found. Please provide a correct file path");
            XDocument doc = XDocument.Load(configFileName);
            initialize(doc, dataConnection, totalNumberOfSpectra, maxTaskCount, scoredAAMods);
        }

        public AnyPTMProcessor(IDataConection dataConnection, CancellationTokenSource m_searchCancel, XDocument doc,
            int RandomSeed, bool UseDiagnosticIons, int maxPeakDepth, int maxIsoformCount, int maxPTMCount,
            ThreeStateEnum scoreNLPeaksToo, double masstolerance, List<AminoAcidModification> scoredAAMods,
            int totalNumberOfSpectra, int maxTaskCount)
        {
            progressMessageQueue = dataConnection.GetProgressMessageQueue();
            if (progressMessageQueue == null)
            {
                throw new ArgumentNullException(
                    "ProgressmessageQueue is null! Please provide a Blocking Collection through the function GetProgressMessageQueue.");
            }
            m_MaximalPeakDepth = maxPeakDepth;
            m_RandomSeed = RandomSeed;
            m_UseDiagnosticIons = UseDiagnosticIons;
            _maxIsoformCount = maxIsoformCount;
            _maxPtmCount = maxPTMCount;
            m_ScoreNLPeaksToo = scoreNLPeaksToo;
            m_MassTolerance = masstolerance;
            ErrorHappened = false;
            m_SearchCancel = m_searchCancel;
            PtmResult = new PTMResultClass();
            _modificationManagment = new ModificationManagement();
            initialize(doc, dataConnection, totalNumberOfSpectra, maxTaskCount, scoredAAMods);
        }

        public void StartPTMLocalisation(PTMResultClass.ResultConsumer resultConsumer)
        {
            try
            {
                Task[] array = new Task[_psmSource.MaxTaskCount];
                progressMessageQueue.Add(new progressMessage(0.0, 0));

                for (int i = 0; i < _psmSource.MaxTaskCount; i++)
                {
                    array[i] = new Task(LocalizePTMs, m_SearchCancel.Token);
                    array[i].Start();
                }
                Task task = new Task(delegate { resultConsumer(PtmResult.IsoformGroupList, m_SearchCancel); },
                    m_SearchCancel.Token, TaskCreationOptions.LongRunning);
                task.Start();
                List<Exception> list = new List<Exception>();
                try
                {
                    Task.WaitAll(array);
                }
                catch (Exception item2)
                {
                    list.Add(item2);
                }
                try
                {
                    PtmResult.CompleteAdding();
                    task.Wait(m_SearchCancel.Token);
                }
                catch (Exception item3)
                {
                    list.Add(item3);
                }
                if (list.Count > 0)
                {
                    bool AllTaskCanceled = false;
                    StringBuilder message = new StringBuilder();
                    foreach (Exception current in list)
                    {
                        if (!(current.GetType() == typeof(OperationCanceledException)) &&
                            current.GetType() == typeof(AggregateException))
                        {
                            AggregateException ex = (AggregateException) current;
                            ex.Flatten()
                                .Handle(delegate(Exception w)
                                {
                                    ErrorHappened = true;
                                    if (w.GetType() == typeof(OperationCanceledException))
                                        return true;
                                    AllTaskCanceled = false;
                                    message.AppendLine(string.Format("{0},\n {1}", w.Message, w.StackTrace));
                                    message.AppendLine("------------\n");
                                    return true;
                                });
                        }
                    }
                    if (AllTaskCanceled)
                    {
                        throw new OperationCanceledException();
                    }
                    progressMessageQueue.TryAdd(new progressMessage(message.ToString()));
                }
            }
            catch (OperationCanceledException)
            {
                throw;
            }
            catch (Exception ex2)
            {
                ErrorHappened = true;
                progressMessageQueue.TryAdd(new progressMessage(ex2.Message + "\n" + ex2.StackTrace));
            }
            finally
            {
                progressMessageQueue.CompleteAdding();
                if (_psmSource != null)
                {
                    _psmSource.Dispose();
                    _psmSource = null;
                }
            }
        }

        private bool comparer(string Sequence, Dictionary<int, Tuple<int, AminoAcidModification>> map1,
            Dictionary<int, Tuple<int, AminoAcidModification>> map2, List<AminoAcidModification> scoredModifications)
        {
            bool result;
            List<char> list = scoredModifications
                .SelectMany(aa => from w in aa.TargetAminoAcids
                    select w.OneLetterCode)
                .ToList();
            foreach (KeyValuePair<int, Tuple<int, AminoAcidModification>> current in map1)
            {
                int key = current.Key;
                int item = current.Value.Item1;
                AminoAcidModification item2 = current.Value.Item2;
                if (scoredModifications.Contains(item2))
                {
                    if (!list.Contains(Sequence[key]))
                    {
                        result = false;
                        return result;
                    }
                    Tuple<int, AminoAcidModification> tuple;
                    if (!map2.TryGetValue(key, out tuple))
                    {
                        result = false;
                        return result;
                    }
                    if (item != tuple.Item1 || !item2.Equals(tuple.Item2))
                    {
                        result = false;
                        return result;
                    }
                }
            }
            foreach (KeyValuePair<int, Tuple<int, AminoAcidModification>> current2 in map2)
            {
                int key2 = current2.Key;
                int item3 = current2.Value.Item1;
                AminoAcidModification item4 = current2.Value.Item2;
                if (scoredModifications.Contains(item4))
                {
                    if (!list.Contains(Sequence[key2]))
                    {
                        result = false;
                        return result;
                    }
                    Tuple<int, AminoAcidModification> tuple;
                    if (!map1.TryGetValue(key2, out tuple))
                    {
                        result = false;
                        return result;
                    }
                    if (item3 != tuple.Item1 || !item4.Equals(tuple.Item2))
                    {
                        result = false;
                        return result;
                    }
                }
            }
            result = true;
            return result;
        }

        private void LocalizePTMs()
        {
            try
            {
                List<SpectraPackageItem> next;
                do
                {
                    next = _psmSource.GetNext();
                    if (next == null)
                        break;

                    if (next.Count != 0)
                    {
                        foreach (SpectraPackageItem current in next)
                        {
                            m_SearchCancel.Token.ThrowIfCancellationRequested();
                            if (current.Error)
                            {
                                PTMResultClass.IsoformGroup isoformGroup =
                                    new PTMResultClass.IsoformGroup(current.ErrorMsg, current.ErrorPeptideIDs);
                                isoformGroup.SequenceID = current.ErrorSequenceID;
                                isoformGroup.SpectrumID = current.ErrorSpectrumID;
                                PtmResult.Add(isoformGroup, m_SearchCancel.Token);
                            }
                            else
                            {
                                PeptideSpectrumMatch psm = current.Psm;
                                string sequence = psm.Sequence.OneLetterCodeString();
                                PeptideScoring peptideScoring = new PeptideScoring(psm, m_MassTolerance,
                                    m_RandomSeed, m_UseDiagnosticIons, 2, m_MaximalPeakDepth,
                                    EquivalentModificationMap, UseMassAccuracyCorrection,
                                    _maxIsoformCount, m_ScoreNLPeaksToo,
                                    getNLsFragmentIonComposition(psm.SpectrumType),
                                    GetFragmentIonTypeComposition(psm.SpectrumType));
                                var isoformListIdPtm = current.IsoformListIdPtm;

                                switch (peptideScoring.Error)
                                {
                                    case PeptideScoring.ScoringError.None:
                                    {
                                        var dictionary =
                                            new Dictionary<PTMSequenceProbability,
                                                Dictionary<int, Tuple<int, AminoAcidModification>>>();
                                        foreach (PTMSequenceProbability current2 in peptideScoring
                                            .PtmSequenceProbabilities)
                                        {
                                            Dictionary<int, Tuple<int, AminoAcidModification>> positionsOf =
                                                current2.Sequence.GetPositionsOf(
                                                    EquivalentModificationMap.ScoredModifications, true);
                                            dictionary.Add(current2, positionsOf);
                                        }
                                        PTMResultClass.IsoformGroup isoformGroup =
                                            new PTMResultClass.IsoformGroup(psm.SequenceOneLetterCode,
                                                peptideScoring.get_SiteProbabilitiesProbabilities(),
                                                peptideScoring.get_Best_SiteProbabilitiesProbabilities(),
                                                peptideScoring.Report);
                                        foreach (SpectraPackageItem.PtmModificationPosition current3 in isoformListIdPtm
                                        )
                                        {
                                            bool flag = false;
                                            foreach (var current4 in dictionary)
                                            {
                                                if (comparer(sequence, current4.Value, current3.ModificationMap,
                                                    EquivalentModificationMap.ScoredModifications))
                                                {
                                                    isoformGroup.add(new PTMResultClass.Peptide(current3.PeptideID,
                                                        current4.Key.Score, current4.Key.Probability));
                                                    flag = true;
                                                    break;
                                                }
                                            }
                                            if (!flag)
                                                isoformGroup.add(new PTMResultClass.Peptide(current3.PeptideID, 0.0,
                                                    0.0));
                                        }
                                        isoformGroup.SequenceID = psm.Sequence.ID;
                                        isoformGroup.SpectrumID = psm.Id;
                                        if (isoformGroup.Count > 0.0)
                                            PtmResult.Add(isoformGroup, m_SearchCancel.Token);
                                        break;
                                    }
                                    case PeptideScoring.ScoringError.TooManyIsoforms:
                                        if (isoformListIdPtm != null)
                                        {
                                            List<int> peptideIDs =
                                                (from w in isoformListIdPtm select w.PeptideID).ToList();
                                            PTMResultClass.IsoformGroup isoformGroup =
                                                new PTMResultClass.IsoformGroup("Too many isoforms", peptideIDs);

                                            isoformGroup.SequenceID = psm.Sequence.ID;
                                            isoformGroup.SpectrumID = psm.Id;
                                            if (isoformGroup.Count > 0.0)
                                                PtmResult.Add(isoformGroup, m_SearchCancel.Token);
                                        }
                                        break;
                                }
                            }
                        }
                    }
                } while (next != null);
            }
            catch (OperationCanceledException)
            {
                throw;
            }
            catch (Exception)
            {
                ErrorHappened = true;
                m_SearchCancel.Cancel();
                throw;
            }
        }

        // 4-24
        public void Dispose()
        {
            Dispose(true);
            GC.SuppressFinalize(this);
        }

        // 4-24
        protected virtual void Dispose(bool disposing)
        {
            if (!disposed)
            {
                if (disposing)
                {
                    if (_psmSource != null)
                    {
                        _psmSource.Dispose();
                        _psmSource = null;
                    }
                    if (m_SearchCancel != null)
                    {
                        m_SearchCancel.Dispose();
                        m_SearchCancel = null;
                    }
                    if (PtmResult != null)
                    {
                        PtmResult.Dispose();
                        PtmResult = null;
                    }
                }
                disposed = true;
            }
        }
    }
}