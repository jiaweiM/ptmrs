using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using ptmrs.Data;

namespace ptmrs
{
    internal class PeptideScoring
    {
        public enum ScoringError
        {
            None,
            NoModificationContained,
            TooManyIsoforms
        }

        private Dictionary<AminoAcidModification, List<Tuple<AminoAcid, int>>> _MaxModNumber;
        private Dictionary<AminoAcidModification, int> _actualModsCount;
        private Dictionary<AminoAcidModification, List<ModificationPosition>> _potentialModPositions;
        private PeakExtractor _peakExtractor;
        private PeptideSpectrumMatch _scoredSpec;
        private List<int> _optimalPeakDepths;
        private List<PTMPepScore> _ptmPepScores = new List<PTMPepScore>();
        private List<AminoAcidModification> _scoredAAMods;
        private ThreeStateEnum _considerNlPeaks;
        private double _massTolerance;
        private double _accuracyCorrection;
        private PTMSiteProbability m_ptmSiteProbs;
        private int _maxPeakDepth = 8;
        private int _minPeakDepth = 2;
        private bool _UseDiagnosticIons;
        private ScoringError _error;
        private StringBuilder _Report;
        private FitComposition _NLsFITComposition;
        private FitComposition _FITComposition;

        private static double _MagicNumber = -10.0 * Math.Log10(0.8);

        public string Report
        {
            get
            {
                List<Tuple<PTMSequenceProbability, List<FragmentIon>>> fIs = GetFIs(true);
                if (fIs.Count >= 2)
                {
                    if (_Report == null)
                    {
                        _Report = new StringBuilder();
                    }
                    else
                    {
                        _Report.Append("; ");
                    }
                    _Report.Append(string.Join("; ", (from w in fIs[0].Item2.Except(fIs[1].Item2)
                        select w.ToString()).ToArray()));
                }
                return _Report?.ToString();
            }
        }

        public ScoringError Error => _error;

        private PeptideSpectrumMatch ScoredSpectrum => _scoredSpec;

        private bool ConsiderNLPeaks => (_scoredSpec.SpectrumType == SpectrumType.HCD ||
                                         _scoredSpec.SpectrumType == SpectrumType.EThcD ||
                                         _considerNlPeaks == ThreeStateEnum.True) &&
                                        _considerNlPeaks != ThreeStateEnum.False;

        public PTMSiteProbability PTMSiteProbabilities => m_ptmSiteProbs;

        /// <summary>
        /// All peptide probabilities. 4-24
        /// </summary>
        public List<PTMSequenceProbability> PtmSequenceProbabilities
        {
            get
            {
                double recPSum = PTMPepScore.GetReciprocalPSum(_ptmPepScores);
                if (recPSum == 0.0)
                    return new List<PTMSequenceProbability>();
                var list = _ptmPepScores.AsParallel()
                    .Select(score => new PTMSequenceProbability(score.Sequence, score.Score,
                        score.ReciprocalScore / recPSum))
                    .ToList();
                list.Sort();
                return list;
            }
        }

        public PeptideScoring(PeptideSpectrumMatch scoredSpec, double massTolerance, int randomSeed,
            bool useDiagnosticIons, int minPeakDepth, int maxPeakDepth, ModificationManagement modificationManagement,
            bool useMassAccuracyCorrection, int maximumIsoformsCount, ThreeStateEnum considerNLPeaks,
            FitComposition nLsFragmentIonComposition, FitComposition fragmentIonComposition)
        {
            _MaxModNumber = modificationManagement.m_maxModificationNumber;
            _scoredSpec = scoredSpec;
            _massTolerance = massTolerance;
            _considerNlPeaks = considerNLPeaks;
            if (useMassAccuracyCorrection)
                _accuracyCorrection = Math.Min(1.0, _massTolerance * 4.0 * 2.0 * 2.0);
            else
                _accuracyCorrection = 1.0;

            _minPeakDepth = minPeakDepth;
            _maxPeakDepth = maxPeakDepth;
            _UseDiagnosticIons = useDiagnosticIons;
            _NLsFITComposition = nLsFragmentIonComposition;
            _FITComposition = fragmentIonComposition;
            _scoredSpec._baseTypes = _scoredSpec.GetBasisFITs(_FITComposition);

            List<AminoAcidModification> scoredModifications = modificationManagement.ScoredModifications;
            _scoredAAMods = new List<AminoAcidModification>(scoredModifications.Count);
            _actualModsCount = new Dictionary<AminoAcidModification, int>(scoredModifications.Count);
            for (int i = 0; i < scoredModifications.Count; i++)
            {
                int numberSites = scoredSpec.Sequence.GetNumberSites(scoredModifications[i]);
                if (numberSites > 0)
                {
                    _actualModsCount.Add(scoredModifications[i], numberSites);
                    _scoredAAMods.Add(scoredModifications[i]);
                }
            }
            _scoredAAMods.TrimExcess();
            if (_actualModsCount.Values.Sum() == 0)
                _error = ScoringError.NoModificationContained;
            else
            {
                double num = 0.0;
                var dictionary = new Dictionary<AminoAcidModification, List<int>>(_scoredAAMods.Count);
                double num2 = 1.0;
                foreach (var current in _scoredAAMods)
                {
                    List<int> positionsOf = _scoredSpec.Sequence.GetPositionsOf(current.TargetAminoAcids,
                        Constants.NO_MOD_ON_C_TERM_K_OR_R, true);
                    dictionary.Add(current, positionsOf);
                    num = Math.Max(PtmMathHelper.BinomialCoefficient(positionsOf.Count, _actualModsCount[current]),num);
                    if (num > maximumIsoformsCount && maximumIsoformsCount != 0)
                    {
                        _error = ScoringError.TooManyIsoforms;
                        return;
                    }
                    if (num > 0.0)
                        num2 *= num;
                    if (num2 > maximumIsoformsCount && maximumIsoformsCount != 0)
                    {
                        _error = ScoringError.TooManyIsoforms;
                        return;
                    }
                }
                _potentialModPositions =
                    new Dictionary<AminoAcidModification, List<ModificationPosition>>();
                List<AminoAcid> unmodifiedAASequence = _scoredSpec.Sequence.GetUnmodifiedAASequence(_scoredAAMods);
                var dictionary2 = new Dictionary<AminoAcidModification, List<List<ModificationPosition>>>();
                List<List<ModificationPosition>> list = new List<List<ModificationPosition>>();
                using (Dictionary<AminoAcidModification, List<int>>.Enumerator enumerator2 = dictionary.GetEnumerator())
                {
                    while (enumerator2.MoveNext())
                    {
                        KeyValuePair<AminoAcidModification, List<int>> modPos = enumerator2.Current;
                        Dictionary<AminoAcidModification, int> arg_2Cc0 = _actualModsCount;
                        KeyValuePair<AminoAcidModification, List<int>> modPos12 = modPos;
                        int num3;
                        arg_2Cc0.TryGetValue(modPos12.Key, out num3);
                        Dictionary<AminoAcidModification, List<Tuple<AminoAcid, int>>> arg_2E80 =
                            _MaxModNumber;
                        KeyValuePair<AminoAcidModification, List<int>> modPos2 = modPos;
                        if (arg_2E80[modPos2.Key].Any(w => w.Item2 > 1))
                        {
                            KeyValuePair<AminoAcidModification, List<int>> modPos3 = modPos;
                            int[] array = new int[modPos3.Value.Count];
                            int num4 = 0;
                            int index;
                            while (true)
                            {
                                int arg_3C20 = num4;
                                KeyValuePair<AminoAcidModification, List<int>> modPos4 = modPos;
                                if (arg_3C20 >= modPos4.Value.Count)
                                    break;
                                KeyValuePair<AminoAcidModification, List<int>> modPos5 = modPos;
                                index = modPos5.Value[num4];
                                AminoAcid aa = unmodifiedAASequence[index];
                                int[] arg_3A40 = array;
                                int arg_3A41 = num4;
                                Dictionary<AminoAcidModification, List<Tuple<AminoAcid, int>>> arg3880 = _MaxModNumber;
                                KeyValuePair<AminoAcidModification, List<int>> modPos6 = modPos;
                                arg_3A40[arg_3A41] = arg3880[modPos6.Key].First(w => w.Item1.Equals(aa)).Item2;
                                num4++;
                            }
                            int arg_3E40 = num3;
                            KeyValuePair<AminoAcidModification, List<int>> modPos7 = modPos;
                            int[][] multiCombinations = PtmMathHelper.GetMultiCombinations(arg_3E40,
                                modPos7.Value.Count, array, maximumIsoformsCount * 4);
                            if (multiCombinations == null || multiCombinations.Length <= 0)
                            {
                                _error = ScoringError.TooManyIsoforms;
                                return;
                            }
                            list = (from possComb in multiCombinations
                                select (from w in possComb.Select((numberOfmodifcations, index) => new
                                    {
                                        numberOfmodifcations,
                                        index
                                    })
                                    where w.numberOfmodifcations > 0
                                    select w).Select(delegate(w)
                                {
                                    KeyValuePair<AminoAcidModification, List<int>> modPos11 = modPos;
                                    return new ModificationPosition(modPos11.Value[w.index], w.numberOfmodifcations,
                                        true);
                                })
                                .ToList<ModificationPosition>()).ToList();
                        }
                        else
                        {
                            KeyValuePair<AminoAcidModification, List<int>> modPos8 = modPos;
                            list = (from w in modPos8.Value.GetCombinationsWithoutReplication(num3)
                                    select (from v in w
                                        select new ModificationPosition(v, 1, true)).ToList())
                                .ToList();
                        }
                        Dictionary<AminoAcidModification, List<List<ModificationPosition>>> arg47B0 = dictionary2;
                        KeyValuePair<AminoAcidModification, List<int>> modPos9 = modPos;
                        arg47B0.Add(modPos9.Key, list);
                        Dictionary<AminoAcidModification, List<ModificationPosition>> arg_4E60 =
                            _potentialModPositions;
                        KeyValuePair<AminoAcidModification, List<int>> modPos10 = modPos;
                        arg_4E60.Add(modPos10.Key, (from w in list.SelectMany(w => w).Distinct<ModificationPosition>()
                            orderby w.ZeroBasedPosition
                            select w).ToList());
                    }
                }
                var possibleCombination = new PtmMathExtensioEnumerableExtensions.Combinations(dictionary2).possibleCombination;
                if (possibleCombination.Count > maximumIsoformsCount && maximumIsoformsCount != 0)
                    _error = ScoringError.TooManyIsoforms;
                else
                {
                    List<AminoAcidSequence> list2 = new List<AminoAcidSequence>(possibleCombination.Count);
                    int num5 = 0;
                    foreach (
                        Dictionary<int, Tuple<int, AminoAcidModification>> current2 in from aaSeq in possibleCombination
                        select aaSeq.SelectMany(aamod =>
                                    from posNumb in aamod.Item2
                                    select new
                                    {
                                        pos = posNumb.OneBasedPosition,
                                        aamodNumb =
                                        new Tuple<int, AminoAcidModification>(posNumb.NumberOfModifications, aamod.Item1)
                                    })
                            .ToDictionary(w => w.pos, w => w.aamodNumb))
                    {
                        num5++;
                        List<AminoAcid> modifiedAASequence = AminoAcidSequence.GetModifiedAASequence(unmodifiedAASequence, current2);
                        AminoAcidSequence item = new AminoAcidSequence(num5, -1, modifiedAASequence,
                            _scoredSpec.Sequence.NTerminus, _scoredSpec.Sequence.CTerminus);
                        list2.Add(item);
                    }
                    list2.TrimExcess();
                    if (randomSeed != -2)
                    {
                        Random random = randomSeed == -1 ? new Random() : new Random(randomSeed);
                        int j = list2.Count;
                        while (j > 1)
                        {
                            j--;
                            int index2 = random.Next(j + 1);
                            AminoAcidSequence value = list2[index2];
                            list2[index2] = list2[j];
                            list2[j] = value;
                        }
                    }
                    UpdatePtmPepScores(list2);
                    UpdatePTMSiteProbs();
                }
            }
        }

        public List<int> GetOptimalPeakDepths()
        {
            return _optimalPeakDepths;
        }

        public List<Tuple<PTMSequenceProbability, List<Tuple<FragmentIon, double>>>> GetFIs(bool matchingOnly)
        {
            List<PTMSequenceProbability> sequenceProbabilities = PtmSequenceProbabilities;
            var list =
                new List<Tuple<PTMSequenceProbability, List<Tuple<FragmentIon, double>>>>(sequenceProbabilities.Count);
            foreach (PTMSequenceProbability probability in sequenceProbabilities)
            {
                AminoAcidSequence sequence = probability.Sequence;
                List<FragmentIonType> allFiTs = PeptideSpectrumMatch.GetAllFiTs(sequence, _scoredAAMods,
                    _scoredSpec._baseTypes, _NLsFITComposition, true);

                List<Tuple<FragmentIon, double>> item = matchingOnly
                    ? PeptideSpectrumMatch.GetMatchingFragmentIons(sequence,
                        _peakExtractor.GetMostIntensePeaks(_optimalPeakDepths), _massTolerance,
                        _scoredSpec.SpectrumType, allFiTs)
                    : sequence.CalcTheoreticalFragIons(_scoredSpec.SpectrumType, allFiTs)
                        .Select(w => new Tuple<FragmentIon, double>(w, -1.0))
                        .ToList();

                list.Add(new Tuple<PTMSequenceProbability, List<Tuple<FragmentIon, double>>>(probability, item));
            }
            return list;
        }

        /// <summary>
        /// 计算肽段打分。4-21
        /// </summary>
        /// <param name="sequence">肽段序列</param>
        /// <param name="peakdepth">当前的 peak depth</param>
        /// <param name="window">窗口索引</param>
        /// <param name="mzWindowBounds">窗口范围</param>
        /// <param name="accuracyCorrection">窗口参数校正因子，默认为1</param>
        /// <returns>肽段打分</returns>
        protected double CalcPeptideScore(AminoAcidSequence sequence, int peakdepth, int window,
            Tuple<double, double> mzWindowBounds, double accuracyCorrection = 1.0)
        {
            AminoAcidSequence.MatchData matchData = sequence.getMatchData(window, peakdepth, _minPeakDepth);
            double probabilityRandamMatch = 2.0 * _massTolerance * matchData.ExtractedPeakCount /
                                            ((mzWindowBounds.Item2 - mzWindowBounds.Item1) * accuracyCorrection);
            int numberTheoreticalFragments;
            switch (_scoredSpec.SpectrumType)
            {
                case SpectrumType.CID_CAD:
                    numberTheoreticalFragments = 8 * (ConsiderNLPeaks ? 2 : 1);
                    break;
                case SpectrumType.ECD_ETD:
                    numberTheoreticalFragments = 6;
                    break;
                case SpectrumType.HCD:
                    numberTheoreticalFragments = 8 * (ConsiderNLPeaks ? 2 : 1);
                    break;
                case SpectrumType.EThcD:
                    numberTheoreticalFragments = 4 + 8 * (ConsiderNLPeaks ? 2 : 1);
                    break;
                default:
                    numberTheoreticalFragments = 20;
                    break;
            }
            return PtmMathHelper.CalculateBinomialScore(probabilityRandamMatch, numberTheoreticalFragments,
                matchData.MatchedPeakCount);
        }

        /// <summary>
        /// 4-21
        /// </summary>
        protected double CalcPeptideScore(AminoAcidSequence sequence, List<int> peakdepth, double mzStart, double mzEnd,
            int consideredFIs, double accuracyCorrection = 1.0)
        {
            AminoAcidSequence.MatchData matchData = sequence.getMatchData(peakdepth, _minPeakDepth);
            double probabilityRandamMatch = 2.0 * _massTolerance * matchData.ExtractedPeakCount / ((mzEnd - mzStart) * accuracyCorrection);
            return PtmMathHelper.CalculateBinomialScore(probabilityRandamMatch, consideredFIs,
                matchData.MatchedPeakCount);
        }

        private void UpdatePtmPepScores(List<AminoAcidSequence> potentialSequences)
        {
            _scoredSpec.ClearPeakAnnotations();
            _peakExtractor = new PeakExtractor(_scoredSpec, _minPeakDepth, _maxPeakDepth);

            if (_UseDiagnosticIons)
            {
                foreach (AminoAcidModification scoredAaMod in _scoredAAMods)
                {
                    if (scoredAaMod?.DiagnosticIons != null)
                    {
                        int numberOfAaModOnSequence = _actualModsCount[scoredAaMod];

                        foreach (DiagnosticIon current2 in from w in scoredAaMod.DiagnosticIons
                            where w.SpectrumType == SpectrumType.None || w.SpectrumType == _scoredSpec.SpectrumType
                            select w)
                        {
                            string value;
                            potentialSequences = current2.FilterPotentialIsoforms(potentialSequences,
                                numberOfAaModOnSequence, scoredAaMod, _scoredSpec.Charge, _massTolerance,
                                _scoredSpec.PrecMassZ, _peakExtractor, out value);
                            if (!string.IsNullOrEmpty(value))
                            {
                                if (_Report == null)
                                    _Report = new StringBuilder();
                                else
                                    _Report.Append("; ");
                                _Report.Append(value);
                            }
                        }
                    }
                }
            }

            if (potentialSequences.Count != 0)
            {
                List<Peak> mostIntensePeaks = _peakExtractor.GetMostIntensePeaks(_maxPeakDepth);

                int num = 0;
                foreach (var aaSeq in potentialSequences)
                {
                    aaSeq.ClearMatchData();
                    List<FragmentIonType> allFiTs = PeptideSpectrumMatch.GetAllFiTs(aaSeq, _scoredAAMods,
                        _scoredSpec._baseTypes, _NLsFITComposition, true);
                    PeptideSpectrumMatch.MatchIsoformToPeaks(aaSeq, mostIntensePeaks, _massTolerance,
                        _scoredSpec.SpectrumType, allFiTs);
                    PeptideSpectrumMatch.InitializeMatchPeaksCount(aaSeq, _peakExtractor, _minPeakDepth,
                        _maxPeakDepth);
                    num = Math.Max(num, aaSeq.ConsiredTheoreticalPeakCount);
                }
                _optimalPeakDepths = new List<int>(_peakExtractor.MzWindowCount);
                List<int> list = Enumerable.Repeat(8, _peakExtractor.MzWindowCount).ToList();
                List<Tuple<int, List<int>>> list2 = new List<Tuple<int, List<int>>>(_peakExtractor.MzWindowCount);

                int num2;
                for (int i = 0; i < _peakExtractor.MzWindowCount; i++)
                {
                    Tuple<double, double> mzWindowBounds = _peakExtractor.GetMzWindowBounds(i);
                    List<PTMPepScoreResults> list3 = new List<PTMPepScoreResults>(potentialSequences.Count);
                    for (int j = 0; j < potentialSequences.Count; j++)
                    {
                        PTMPepScoreResults pTmPepScoreResults = new PTMPepScoreResults(potentialSequences[j]);
                        for (int k = _minPeakDepth; k <= _maxPeakDepth; k++)
                        {
                            double score = CalcPeptideScore(potentialSequences[j], k, i, mzWindowBounds,
                                _accuracyCorrection);
                            pTmPepScoreResults.AddScoreResult(k, score);
                        }
                        list3.Add(pTmPepScoreResults);
                    }
                    List<int> item;
                    num2 = FindOptimalPeakDepth(list3, out item);
                    list2.Add(Tuple.Create(i, item));
                    list[i] = num2;
                    _optimalPeakDepths.Add(num2);
                }

                var tuples = list2.Where(depth => depth.Item2.Count >= 2);
                foreach (var tuple in tuples)
                {
                    List<PTMPepScoreResults> list3 = new List<PTMPepScoreResults>(potentialSequences.Count);
                    var num3 = -1.0;
                    var value2 = -1;
                    for (int l = 0; l < potentialSequences.Count; l++)
                    {
                        PTMPepScoreResults pTmPepScoreResults = new PTMPepScoreResults(potentialSequences[l]);
                        foreach (int current5 in tuple.Item2)
                        {
                            list[tuple.Item1] = current5;
                            double num4 = CalcPeptideScore(potentialSequences[l], list, _peakExtractor.StartMassZ,
                                _peakExtractor.EndMassZ, num, _accuracyCorrection);
                            pTmPepScoreResults.AddScoreResult(current5, num4);
                            if (num4 > num3)
                            {
                                num3 = num4;
                                value2 = current5;
                            }
                        }
                        list3.Add(pTmPepScoreResults);
                    }
                    list[tuple.Item1] = value2;
                    _optimalPeakDepths[tuple.Item1] = value2;
                }

                _ptmPepScores = new List<PTMPepScore>(potentialSequences.Count);
                foreach (AminoAcidSequence aaSeq in potentialSequences)
                {
                    double peptideScore = CalcPeptideScore(aaSeq, _optimalPeakDepths,
                        _peakExtractor.StartMassZ, _peakExtractor.EndMassZ, num, _accuracyCorrection);
                    _ptmPepScores.Add(new PTMPepScore(aaSeq, peptideScore));
                }
            }
        }

        /// <summary>
        /// 查找当前窗口最佳的 peakdepth 4.20
        /// </summary>
        /// <param name="wndPepScores">当前窗口所有的打分值</param>
        /// <param name="outEquivalentPeakDepths">输出的所有最佳 peak depth</param>
        /// <returns>最佳 peak depth</returns>
        private int FindOptimalPeakDepth(List<PTMPepScoreResults> wndPepScores, out List<int> outEquivalentPeakDepths)
        {
            outEquivalentPeakDepths = new List<int>();

            // 添加所有可能的 peak depth
            List<int> equivalentPeakDepths = new List<int>(_maxPeakDepth - _minPeakDepth + 1);
            for (int i = _minPeakDepth; i <= _maxPeakDepth; i++)
                equivalentPeakDepths.Add(i);

            // 如果当前窗口只有一个肽段的打分值，则输出最高打分值处的 peak depth
            if (wndPepScores.Count <= 1)
            {
                outEquivalentPeakDepths.AddRange(equivalentPeakDepths);
                int depth;
                wndPepScores[0].GetMaxScoreResult(_minPeakDepth, out depth);
                return depth;
            }

            int peakDepth = _minPeakDepth;
            int itCount = wndPepScores.Count <= 3 ? wndPepScores.Count - 1 : 3;
            for (int it = 1; it <= itCount; it++)
            {
                peakDepth = optimizePeakDepth(wndPepScores, it, ref equivalentPeakDepths);
                if (equivalentPeakDepths.Count <= 1)
                    break;
            }
            outEquivalentPeakDepths.AddRange(equivalentPeakDepths);
            return peakDepth;
        }

        /// <summary>
        /// 查找当前 iteration 下的最佳 peak depth 4.20
        /// </summary>
        /// <param name="wndPepScores"></param>
        /// <param name="iteration">对比的次序，如果是1，则对比排序第1和第2肽段</param>
        /// <param name="equivalentPeakDepths">所有的 peak depth，从其中找出最佳值</param>
        /// <returns>best depth</returns>
        private int optimizePeakDepth(List<PTMPepScoreResults> wndPepScores, int iteration,
            ref List<int> equivalentPeakDepths)
        {
            int bestDepth = -1;
            double maxDeltaScore = -1.0;
            double maxScore = -1.0;
            List<int> list = new List<int>(equivalentPeakDepths.Count);
            foreach (int depth in equivalentPeakDepths)
            {
                PTMPepScoreResults.SortPTMPeptideScores(wndPepScores, depth);
                double deltaScore = Math.Abs(wndPepScores[iteration - 1].GetScoreResult(depth) -
                                             wndPepScores[iteration].GetScoreResult(depth));
                if (deltaScore > maxDeltaScore)
                {
                    maxDeltaScore = deltaScore;
                    if (list.Count != 0)
                        list.Clear();
                    list.Add(depth);
                    maxScore = wndPepScores[iteration - 1].GetScoreResult(depth);
                    bestDepth = depth;
                }
                else if (PtmMathHelper.Equal(deltaScore, maxDeltaScore, 0.001))
                {
                    list.Add(depth);
                    if (wndPepScores[iteration - 1].GetScoreResult(depth) > maxScore)
                    {
                        maxScore = wndPepScores[iteration - 1].GetScoreResult(depth);
                        bestDepth = depth;
                    }
                }
            }
            equivalentPeakDepths = list;
            return bestDepth;
        }

        protected void UpdatePTMSiteProbs()
        {
            double reciprocalPSum = PTMPepScore.GetReciprocalPSum(_ptmPepScores);
            if (Math.Abs(reciprocalPSum) > 0.001)
            {
                double num;
                switch (Constants.SITEPROBABILITY_MODE)
                {
                    case PTMSiteProbMode.AtLeastNTime:
                    {
                        List<PTMSiteProbability.SiteProbability> list =
                            new List<PTMSiteProbability.SiteProbability>(
                                _potentialModPositions.Sum(w => w.Value.Max(v => v.NumberOfModifications)));
                        using (
                            Dictionary<AminoAcidModification, List<ModificationPosition>>.Enumerator enumerator =
                                _potentialModPositions.GetEnumerator())
                        {
                            while (enumerator.MoveNext())
                            {
                                KeyValuePair<AminoAcidModification, List<ModificationPosition>> potPositions =
                                    enumerator.Current;
                                KeyValuePair<AminoAcidModification, List<ModificationPosition>> potPositions5 =
                                    potPositions;
                                potPositions5.Value.Max(w => w.NumberOfModifications);
                                KeyValuePair<AminoAcidModification, List<ModificationPosition>> potPositions2 =
                                    potPositions;
                                foreach (IGrouping<int, ModificationPosition> current in from w in potPositions2.Value
                                    group w by w.ZeroBasedPosition)
                                {
                                    int pos = current.Key;
                                    foreach (ModificationPosition current2 in current)
                                    {
                                        int numberOfModification = current2.NumberOfModifications;
                                        num = _ptmPepScores.Where(delegate(PTMPepScore pepScore)
                                            {
                                                AminoAcidSequence arg_2A0 = pepScore.Sequence;
                                                int arg_2A1 = pos;
                                                int arg_2A2 = numberOfModification;
                                                KeyValuePair<AminoAcidModification, List<ModificationPosition>>
                                                    potPositions4 = potPositions;
                                                return arg_2A0.HasModification(arg_2A1, arg_2A2, potPositions4.Key);
                                            })
                                            .Sum(w => w.ReciprocalScore);
                                        List<PTMSiteProbability.SiteProbability> arg_1C30 = list;
                                        int arg1Be0 = pos;
                                        int arg1Be1 = numberOfModification;
                                        KeyValuePair<AminoAcidModification, List<ModificationPosition>> potPositions3 =
                                            potPositions;
                                        arg_1C30.Add(new PTMSiteProbability.SiteProbability(arg1Be0, arg1Be1,
                                            potPositions3.Key, num / reciprocalPSum));
                                    }
                                }
                            }
                        }
                        m_ptmSiteProbs = new PTMSiteProbability(list);
                        break;
                    }
                    case PTMSiteProbMode.ExactlyNTime:
                    {
                        List<PTMSiteProbability.SiteProbability> list =
                            new List<PTMSiteProbability.SiteProbability>(
                                _potentialModPositions.SelectMany(w => from v in w.Value
                                        select v.NumberOfModifications)
                                    .Sum());
                        using (
                            Dictionary<AminoAcidModification, List<ModificationPosition>>.Enumerator enumerator4 =
                                _potentialModPositions.GetEnumerator())
                        {
                            while (enumerator4.MoveNext())
                            {
                                KeyValuePair<AminoAcidModification, List<ModificationPosition>> potentialPosition =
                                    enumerator4.Current;
                                KeyValuePair<AminoAcidModification, List<ModificationPosition>> potentialPosition4 =
                                    potentialPosition;
                                using (
                                    List<ModificationPosition>.Enumerator enumerator5 =
                                        potentialPosition4.Value.GetEnumerator())
                                {
                                    while (enumerator5.MoveNext())
                                    {
                                        ModificationPosition positions = enumerator5.Current;
                                        num = 0.0;
                                        foreach (
                                            PTMPepScore current3 in
                                            _ptmPepScores.Where(delegate(PTMPepScore potentialAASeq)
                                            {
                                                AminoAcidSequence arg_1F0 = potentialAASeq.Sequence;
                                                ModificationPosition arg_1F1 = positions;
                                                KeyValuePair<AminoAcidModification, List<ModificationPosition>>
                                                    potentialPosition3 = potentialPosition;
                                                return arg_1F0.HasModification(arg_1F1, potentialPosition3.Key);
                                            }))
                                        {
                                            num += current3.ReciprocalScore;
                                        }
                                        List<PTMSiteProbability.SiteProbability> arg34B0 = list;
                                        ModificationPosition arg3460 = positions;
                                        KeyValuePair<AminoAcidModification, List<ModificationPosition>>
                                            potentialPosition2 = potentialPosition;
                                        arg34B0.Add(new PTMSiteProbability.SiteProbability(arg3460,
                                            potentialPosition2.Key,
                                            num / reciprocalPSum / positions.NumberOfModifications));
                                    }
                                }
                            }
                        }
                        m_ptmSiteProbs = new PTMSiteProbability(list);
                        break;
                    }
                    default:
                        throw new ArgumentException("Unknown site probability count mode");
                }
            }
        }

        public Dictionary<AminoAcidModification, string> get_SiteProbabilitiesProbabilities()
        {
            PTMSiteProbability siteProbabilities = PTMSiteProbabilities;
            return _scoredAAMods.ToDictionary(w => w,
                w => siteProbabilities.GetSiteProbabilityStrings(_scoredSpec.SequenceOneLetterCode, w));
        }

        public string get_Best_SiteProbabilitiesProbabilities()
        {
            string result;
            double num = _ptmPepScores.Sum(w => w.Score);
            if (num == 0.0)
            {
                result = "Inconclusive data";
            }
            else
            {
                _ptmPepScores.Sort();
                PTMSiteProbability pTMSiteProbabilities = PTMSiteProbabilities;
                double maxScore = _ptmPepScores.Max(w => w.Score);
                List<Tuple<AminoAcidModification, List<AminoAcidSequence>>> aminoAcidModificationIsoformlists =
                (from aamod in _scoredAAMods
                    select
                    new Tuple<AminoAcidModification, List<AminoAcidSequence>>(aamod,
                        (from w in
                            _ptmPepScores.GroupBy(w => w.Sequence.GetModificationMap(aamod),
                                new ModificationPosition.Compare())
                            where maxScore - w.Max(v => v.Score) <= _MagicNumber
                            select w
                            into v
                            select v.First().Sequence).ToList())).ToList();
                result = pTMSiteProbabilities.GetSiteProbabilityString(_scoredSpec.SequenceOneLetterCode,
                    aminoAcidModificationIsoformlists);
            }
            return result;
        }
    }
}