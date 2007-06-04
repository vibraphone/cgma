extern "C" int c_main(int argc, char **argv);

int main(int argc, char **argv)
{
  int result = c_main(argc, argv);
  return result;
}
