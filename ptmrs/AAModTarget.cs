using System;

namespace ptmrs
{
    // 4.19
    public abstract class AAModTarget
    {
        protected AminoAcidModification m_mod;
        protected int m_numberOfModifications;

        public bool IsModified => m_mod != null;

        public AminoAcidModification Modification => m_mod;

        public int NumberOfModifications => m_numberOfModifications;

        public abstract double GetTotalMass();

        public abstract void setModification(AminoAcidModification modification, int number);

        public abstract void setModification(Tuple<int, AminoAcidModification> modification);
    }
}