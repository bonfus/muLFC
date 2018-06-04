include(CheckCXXSourceCompiles)

# Taken from: https://github.com/elemental/
#
# C++11 random number generation
# ==============================
# Note: It was noticed that, for certain relatively recent Intel compiler
#       releases, the following snippets will suggest that normal distributions
#       are supported, but not uniform integer or real distributions. 
#       Unfortunately, it has been observed that, when one attempts to sample
#       from a normal distribution, the application hangs indefinitely.
#       Thus, Elemental takes the approach of requiring all of the following
#       to be detected before any are used.
set(NORMAL_CODE
    "#include <random>
     int main()
     {
         std::random_device rd;
         std::mt19937 mt(rd());
         std::normal_distribution<double> dist(0,1);
         const double x = dist(mt);
         return 0;
     }")
set(UNIFORM_INT_CODE
    "#include <random>
     int main()
     {
         std::random_device rd;
         std::mt19937 mt(rd());
         std::uniform_int_distribution<int> dist(0,10);
         const int x = dist(mt);
         return 0;
     }")
set(UNIFORM_REAL_CODE
    "#include <random>
     int main()
     {
         std::random_device rd;
         std::mt19937 mt(rd());
         std::uniform_real_distribution<double> dist(0,1);
         const double x = dist(mt);
         return 0;
     }")
check_cxx_source_compiles("${NORMAL_CODE}" HAVE_NORMAL_DIST)
check_cxx_source_compiles("${UNIFORM_INT_CODE}" HAVE_UNIFORM_INT_DIST)
check_cxx_source_compiles("${UNIFORM_REAL_CODE}" HAVE_UNIFORM_REAL_DIST)
if(HAVE_NORMAL_DIST AND 
   HAVE_UNIFORM_INT_DIST AND
   HAVE_UNIFORM_REAL_DIST)
  set(HAVE_CXX11RANDOM TRUE)
else()
  set(HAVE_CXX11RANDOM FALSE)
endif()
