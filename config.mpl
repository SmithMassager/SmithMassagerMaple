# Number of threads to use for maple
TASKSIZE := 40;

# Enable external c library
if (FileTools[Exists]("lib") and FileTools[IsDirectory]("lib")) then
  CLIB := true;
  print("Using C libraries");
else
  CLIB := false;
  print("NOT using C libraries");
end if;
