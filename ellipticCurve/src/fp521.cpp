#ifdef YKM_ECC_USE_MPN
#include "FP-mpn-impl.hpp"
size_t Fp::size_ = 9;

#else
#include "FP-impl.hpp"

#endif
