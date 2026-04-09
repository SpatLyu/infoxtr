#ifndef INFOXTR_KOCMI_HPP
#define INFOXTR_KOCMI_HPP

#include <vector>
#include <cmath>
#include <limits>
#include <thread>
#include <cstdint>
#include <numeric>
#include <algorithm>
#include <stdexcept>
#include "infoxtr/numericutils.hpp"
#include "infoxtr/distance.hpp"
#include "infoxtr/infotheo.hpp"
#include <RcppThread.h>

namespace infoxtr 
{

namespace kocmi 
{
    using DiscVec = std::vector<uint64_t>;
    using DiscMat = std::vector<std::vector<uint64_t>>;
    using ContVec = std::vector<double>;
    using DiscMat = std::vector<std::vector<double>>;

} // namespace kocmi

}

#endif // INFOXTR_KOCMI_HPP
