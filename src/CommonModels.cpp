/* Author: Masaki Murooka */

#include <CCC/CommonModels.h>
#include <CCC/Constants.h>

using namespace CCC;

ComZmpModelJerkInput::ComZmpModelJerkInput(double com_height)
{
  A_(0, 1) = 1;
  A_(1, 2) = 1;

  B_(2, 0) = 1;

  C_(0, 0) = 1;
  C_(0, 2) = -1 * com_height / constants::g;
}
