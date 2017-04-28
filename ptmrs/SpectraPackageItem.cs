using System;
using System.Collections.Generic;
using System.Linq;

namespace ptmrs
{
    public class SpectraPackageItem
    {
        private readonly double _spectraCount;

        public struct PtmModificationPosition
        {
            public int PeptideID;
            public Dictionary<int, Tuple<int, AminoAcidModification>> ModificationMap;

            public PtmModificationPosition(int peptideId,
                Dictionary<int, Tuple<int, AminoAcidModification>> modificationMap)
            {
                PeptideID = peptideId;
                ModificationMap = modificationMap;
            }

            // 4-24
            public PtmModificationPosition(int peptideID, string modMap, List<AminoAcidModification> modifications)
            {
                PeptideID = peptideID;
                if (string.IsNullOrEmpty(modMap))
                    throw new ArgumentException("Modification map is empty");

                var dictionary = new Dictionary<int, Tuple<int, AminoAcidModification>>();

                var modNames = modMap.Split(";".ToCharArray(), StringSplitOptions.RemoveEmptyEntries).Select(w => w.Split(",".ToCharArray()));
                foreach (var name in modNames)
                {
                    int key;
                    int item;
                    int modId;
                    if (!int.TryParse(name[0], out key) || !int.TryParse(name[1], out item) ||
                        !int.TryParse(name[2], out modId))
                        throw new ArgumentException($"Parse to integer error: {string.Join(" ", name)}");
                    AminoAcidModification aminoAcidModification = modifications.Find(w => w.ID == modId);
                    if (aminoAcidModification == null)
                        throw new ArgumentException("Modification not found (compared by ID)");

                    dictionary.Add(key, new Tuple<int, AminoAcidModification>(item, aminoAcidModification));
                }
                ModificationMap = dictionary;
            }
        }

        public List<PtmModificationPosition> IsoformListIdPtm;

        public PeptideSpectrumMatch Psm;

        private string _errorMsg;

        private List<int> _errorPeptideIds;

        private int _errorSpectrumID;

        private int _errorSequenceID;

        public bool Error { get; private set; }

        public string ErrorMsg
        {
            get
            {
                if (Error)
                {
                    return _errorMsg;
                }
                return "";
            }
            private set => _errorMsg = value;
        }

        public List<int> ErrorPeptideIDs
        {
            get
            {
                if (Error)
                {
                    return _errorPeptideIds;
                }
                return new List<int>();
            }
            private set => _errorPeptideIds = value;
        }

        public int ErrorSpectrumID
        {
            get
            {
                if (Error)
                {
                    return _errorSpectrumID;
                }
                return -1;
            }
            private set => _errorSpectrumID = value;
        }

        public int ErrorSequenceID
        {
            get
            {
                if (Error)
                {
                    return _errorSequenceID;
                }
                return -1;
            }
            set => _errorSequenceID = value;
        }

        public int PeptideCount
        {
            get
            {
                if (_errorPeptideIds != null)
                {
                    return _errorPeptideIds.Count;
                }
                if (IsoformListIdPtm != null)
                {
                    return IsoformListIdPtm.Count;
                }
                return 0;
            }
        }

        public double SpectraCount { get; private set; }

        public SpectraPackageItem(List<int> peptideIds, double spectraCount, string msg, int spectrumID, int sequenceID)
        {
            _spectraCount = spectraCount;
            ErrorMsg = msg;
            Error = true;
            ErrorPeptideIDs = peptideIds;
            ErrorSequenceID = sequenceID;
            ErrorSpectrumID = spectrumID;
        }

        public SpectraPackageItem(PeptideSpectrumMatch psm, double spectraCount,
            List<PtmModificationPosition> isoformListIdPtm)
        {
            try
            {
                Psm = psm ?? throw new ArgumentNullException("PSM is null");
                SpectraCount = spectraCount;
                IsoformListIdPtm = isoformListIdPtm;
                Error = false;
                ErrorMsg = "";
            }
            catch (Exception ex)
            {
                Error = true;
                ErrorMsg = ex.Message;
                throw;
            }
        }
    }
}