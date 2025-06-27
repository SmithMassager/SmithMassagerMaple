# Number of threads to use for maple
TASKSIZE := 40;

# Enable external c library
if (FileTools[Exists]("lib") and FileTools[IsDirectory]("lib")) then
  CLIB := true;
else
  CLIB := false;
end if;

print(CLIB);
