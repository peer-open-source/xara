
constexpr unsigned long log2(unsigned long n)
{
  return ( (n<2) ? 1 : 1+log2(n/2));
}


class FrameSection {
};


#if 0
enum class ShellType {
  N, Vy, Vz,
  T, My, Mz
};

class TrussSection {
};

class ShellSection {
};

#endif
