#ifndef CFF_IO_H
#define CFF_IO_H

#include "generate_cff.h"

int  write_cff(const CFF *cff, const char *filename);
CFF *read_cff(const char *filename);

#endif /* CFF_IO_H */
