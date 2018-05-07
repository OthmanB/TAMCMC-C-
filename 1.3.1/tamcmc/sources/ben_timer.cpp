#include "unistd.h"
#include <sys/time.h>

long ben_clock()
{
	struct timeval tnow;
	gettimeofday(&tnow,0);
	return(tnow.tv_sec);
}
