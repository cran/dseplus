--
-- test function
--
function compare.f5.annual

  set compare.f5.count = compare.f5.count + 1
  if 2*trunc(compare.f5.count/2) ne compare.f5.count
    return (compare.c5.annual/100000 + compare.c5.annual)/100000
  else
    return 100 + (compare.c5.annual/100000 + compare.c5.annual)/100000
  endif

end function
