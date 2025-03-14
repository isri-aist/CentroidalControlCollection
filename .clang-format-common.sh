# This script is meant to be sourced from other scripts

# Check for clang-format, prefer 14 if available
if [[ -x "$(command -v clang-format-14)" ]]; then
  clang_format=clang-format-14
elif [[ -x "$(command -v clang-format)" ]]; then
  clang_format=clang-format
else
  echo "clang-format or clang-format-14 must be installed"
  exit 1
fi

# Find all source files in the project minus those that are auto-generated or we do not maintain
src_files=`find include src tests -name '*.cpp' -or -name '*.h' -or -name '*.hpp'`
