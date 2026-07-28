__config() -> {
    'commands' -> {
        '<chunk_radius>' -> 'geodegen'
    },
    'arguments' -> {
        'chunk_radius' -> {'type' -> 'int', 'min' -> 1, 'max' -> 64, 'suggest' -> [6]}
    }
};

geodegen(chunk_radius) -> (
    block_radius = 16 * (chunk_radius * 2 + 1) / 2;

    // determine y bounds
    game_version = system_info('game_version');
    min_y = -64;
    max_y = 45;

    if(game_version == '1.17' || game_version == '1.17.1',
        min_y = 0;
        max_y = 60;
    );

    finder_output = read_file('output', 'shared_json');
    for(finder_output,
        generation_task(_:'center_x', _:'center_z', min_y, max_y, block_radius);
    );

    print('Finished!');
);


generation_task(x, z, min_y, max_y, radius) -> (
    // geodes can leak out of a chunk about this much
    bonus = radius + 8;
    count = volume((x - bonus, min_y, z - bonus), (x + bonus, max_y, z + bonus),
      _ == 'budding_amethyst'
    );
    print(str('%d budding amethyst centered at (%d, %d)', count, x, z));
);
