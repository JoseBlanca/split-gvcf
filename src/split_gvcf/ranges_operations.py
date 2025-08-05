class RangeError(ValueError):
    pass


class BeforeFirstRange(RangeError):
    pass


class InFirstRange(RangeError):
    pass


def bisect_left(a, x, lo=0, hi=None, *, key=None):
    """Return the index where to insert item x in list a, assuming a is sorted.

    The return value i is such that all e in a[:i] have e < x, and all e in
    a[i:] have e >= x.  So if x already appears in the list, a.insert(i, x) will
    insert just before the leftmost x already there.

    Optional args lo (default 0) and hi (default len(a)) bound the
    slice of a to be searched.

    A custom key function can be supplied to customize the sort order.
    """

    if lo < 0:
        raise ValueError("lo must be non-negative")
    if hi is None:
        hi = len(a)
    # Note, the comparison uses "<" to match the
    # __lt__() logic in list.sort() and in heapq.
    if key is None:
        while lo < hi:
            mid = (lo + hi) // 2
            if a[mid] < x:
                lo = mid + 1
            else:
                hi = mid
    else:
        while lo < hi:
            mid = (lo + hi) // 2
            if key(a[mid]) < x:
                lo = mid + 1
            else:
                hi = mid
    return lo


def _add_seqname(seqname, sorted_seqnames, seqnames, first_range_for_seq, idx):
    if seqname in seqnames:
        raise ValueError(
            "Trying to add a seqname that was already found, maybe they are not in order"
        )
    sorted_seqnames.append(seqname)
    seqnames.add(seqname)
    first_range_for_seq[seqname] = idx


class RangesSearcher:
    def __init__(self, ranges):
        self.ranges = ranges
        self.seqnames = ranges.get_seqnames()

        sorted_seqnames = []
        seqnames = set()
        current_seq = None
        first_range_for_seq = {}
        for idx, seqname in enumerate(self.seqnames):
            if current_seq is None or seqname != current_seq:
                _add_seqname(
                    seqname, sorted_seqnames, seqnames, first_range_for_seq, idx
                )
                current_seq = seqname
        self.sorted_seqnames = sorted_seqnames
        self.seqname_order = {
            seqname: idx for idx, seqname in enumerate(sorted_seqnames)
        }

        self.starts = ranges.get_start()
        self.ends = ranges.get_end()
        self.first_range_for_seq = first_range_for_seq

    def _before_than(self, seqname1, pos1, seqname2, pos2):
        if seqname1 != seqname2:
            return self.seqname_order[seqname1] < self.seqname_order[seqname2]
        return pos1 < pos2

    def find_prev_range(self, seqname, pos):
        lo = 0
        hi = len(self.starts)

        while lo < hi:
            mid = (lo + hi) // 2
            if self._before_than(self.seqnames[mid], self.ends[mid], seqname, pos):
                lo = mid + 1
            else:
                hi = mid

        if lo >= len(self.ranges):
            return len(self.ranges) - 1

        # the position might be inside the found range
        if lo > 0 and self._before_than(seqname, pos, self.seqnames[lo], self.ends[lo]):
            lo -= 1

        first_range_for_seq = self.first_range_for_seq[seqname]
        if seqname != self.seqnames[lo]:
            # We are before or in the first range of a chrom and the range belongs to the prev. chrom
            lo += 1
        if lo == first_range_for_seq:
            if self._before_than(seqname, pos, self.seqnames[lo], self.starts[lo]):
                raise BeforeFirstRange()
            elif self._before_than(seqname, pos, self.seqnames[lo], self.ends[lo] + 1):
                raise InFirstRange()

        return lo
