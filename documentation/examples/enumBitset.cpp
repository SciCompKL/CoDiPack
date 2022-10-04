
#include <iostream>
#include <codi.hpp>

//! [Enum definition]
enum class TestFlags {
  Bit1 = 0,
  Bit2,
  Bit3,
  Bit4,
  Bit5,
  MaxElement,
};

using TestOptions = codi::EnumBitset<TestFlags>;

CODI_INLINE TestOptions operator|(TestFlags a, TestFlags b) {
  return TestOptions(a) | b;
}
//! [Enum definition]



int main(int nargs, char** args) {
//! [Enum use]
  TestOptions options = TestFlags::Bit1 | TestFlags::Bit3 | TestFlags::Bit5;

  std::cout << "Flags: " << options << std::endl;

  std::cout << "Flag1: " << options.test(TestFlags::Bit1) << std::endl;
  std::cout << "Flag2: " << options.test(TestFlags::Bit2) << std::endl;
  std::cout << "Flag3: " << options.test(TestFlags::Bit3) << std::endl;
  std::cout << "Flag4: " << options.test(TestFlags::Bit4) << std::endl;
  std::cout << "Flag5: " << options.test(TestFlags::Bit5) << std::endl;
//! [Enum use]
}
