#include <stdlib.h>
#include <pwd.h>
#include <sys/types.h>
#include <unistd.h>

int main(int argc, char**argv)
{
  struct passwd *ps;
  char    *svcuser;

  ps = getpwuid(geteuid());
  svcuser = ps->pw_name;
  printf("%s\n", svcuser);

  return 0;
}
