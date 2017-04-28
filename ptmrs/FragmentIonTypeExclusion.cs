using System;
using System.Collections.Generic;

namespace ptmrs
{
    public class FragmentIonTypeExclusion
    {
        private int m_curExclusionListIndex;

        private List<Tuple<FragmentIonType, List<int>>> m_exclusionLists =
            new List<Tuple<FragmentIonType, List<int>>>();

        public FragmentIonTypeExclusion()
        {
            m_curExclusionListIndex = -1;
            m_exclusionLists = new List<Tuple<FragmentIonType, List<int>>>();
        }

        public void AddExclusionType(FragmentIonType type)
        {
            AddExclusionList(new Tuple<FragmentIonType, List<int>>(type, new List<int>()));
        }

        public void AddExclusionIndexToCurrentType(int fiIndex)
        {
            m_exclusionLists[m_curExclusionListIndex].Item2.Add(fiIndex);
        }

        public void AddExclusionList(Tuple<FragmentIonType, List<int>> exclusionList)
        {
            m_exclusionLists.Add(exclusionList);
            m_curExclusionListIndex = m_exclusionLists.Count - 1;
        }

        public bool IsExclusionTypeToAdd(FragmentIonType type)
        {
            foreach (Tuple<FragmentIonType, List<int>> t in m_exclusionLists)
            {
                if (t.Item1.Equals(type, true, true, true, true, true))
                {
                    return false;
                }
            }
            return true;
        }

        public bool CheckFragmentIonIndex(int fiIndex)
        {
            List<int> item = m_exclusionLists[m_curExclusionListIndex].Item2;
            foreach (int t in item)
            {
                if (t == fiIndex)
                {
                    return true;
                }
            }
            return false;
        }

        public bool IsEmpty()
        {
            return m_exclusionLists == null || m_exclusionLists.Count == 0;
        }
    }
}