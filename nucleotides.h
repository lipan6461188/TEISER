#define _U				0x01
#define _C				0x02
#define _G				0x04
#define _A				0x08
#define _N				0x0F

#define _Y				0x03 //UC
#define _R				0x0C //AG
#define _K				0x05 //UG
#define _M				0x0A //AC
#define _S				0x06 //GC
#define _W				0x09 //AU

#define _B				0x07 //GUC
#define _D				0x0D //GAU
#define _H				0x0B //ACU
#define _V				0x0E //GCA

#define NUM_LETTERS		4
#define LAST_LETTER		_A

static char lcl_base_templates1[] = "TCGA";
static char lcl_base_templates2[] = "tcga";
static char lcl_base_templates3[] = "UCGA";
static char lcl_base_templates4[] = "ucga";

static char lcl_N1 = 'N';
static char lcl_N2 = 'n';

#define _pair			1
#define _leftBulge		2
#define _rightBulge		3
