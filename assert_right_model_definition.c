#include <MODEL.h>

void assert_right_model_definition( Parameter_Table * P )
{

#if defined  DIFFUSION_HII_nD
  
    assert ( P->TYPE_of_MODEL == 16 );
  
#endif
}
