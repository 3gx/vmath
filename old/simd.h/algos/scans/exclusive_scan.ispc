static int warp_scan(int value)
{
  for (int i = 1; i < programCount; i *= 2)
  {
    int n = shuffle(value, max(programIndex-i,0));
    if (programIndex >= i)
      value += n;
  }
  return value;
}

static int warp_scan16(int value)
{
  int n, i;
  i = 1; value += shuffle(value, programIndex-i) & (-(int)(programIndex >= i));
  i = 2; value += shuffle(value, programIndex-i) & (-(int)(programIndex >= i));
  i = 4; value += shuffle(value, programIndex-i) & (-(int)(programIndex >= i));
  if (programCount >= 8)
  {
    i = 8; value += shuffle(value, programIndex-i) & (-(int)(programIndex >= i));
    if (programCount >= 16) 
      i = 16; value += shuffle(value, programIndex-i) & (-(int)(programIndex >= i));
  }
  return value;
}

export void exclusive_scan(
    uniform const int n,
    uniform       int data[])
{
  int carryValue = 0;
  int x;
  int y = 0;
  foreach(i=0...n)
  {
    carryValue += extract(y,programCount-1);
    x = data[i];
    y = warp_scan16(x);
    data[i] = y + carryValue - x;
  }
}

export void add(
    uniform const int n,
    uniform const int value,
    uniform       int data[])
{
  foreach(i=0...n)
  {
    data[i] += value;
  }
}

