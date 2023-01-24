__command() -> (
    finder_output = read_file('geodes', 'shared_text');
    for(finder_output, 
         x = number(_ ~ 'at ([-]*\\d+)') * 1;
         z = number(_ ~ '([-]*\\d+)$') * 1;
         generation_task(x,z,10);
         print('Generated area surrounding '+x+', '+z);
    );
    print('Finished');
);

generation_task(x, z, rad) -> (
  for(range(-rad, rad+1),
    dx = _;
    for(range(-rad, rad+1), 
      dz = _;
      block(x+dx*16, 0, z+dz*16);
    )
  );
)
