// INFO: Configure as the entrypoint in a Fabric project to verify whether new MC versions will need
// to be added to geode-finder. May need modifications.

package geode;

import it.unimi.dsi.fastutil.doubles.DoubleList;
import java.lang.reflect.Proxy;
import java.util.List;
import java.util.concurrent.atomic.AtomicInteger;
import net.fabricmc.api.ModInitializer;
import net.fabricmc.fabric.api.event.lifecycle.v1.ServerLifecycleEvents;
import net.minecraft.core.BlockPos;
import net.minecraft.core.Registry;
import net.minecraft.resources.ResourceKey;
import net.minecraft.resources.ResourceLocation;
import net.minecraft.server.MinecraftServer;
import net.minecraft.util.Mth;
import net.minecraft.world.level.WorldGenLevel;
import net.minecraft.world.level.biome.BiomeSource.StepFeatureData;
import net.minecraft.world.level.block.Blocks;
import net.minecraft.world.level.block.state.BlockState;
import net.minecraft.world.level.chunk.ChunkGenerator;
import net.minecraft.world.level.levelgen.LegacyRandomSource;
import net.minecraft.world.level.levelgen.RandomSource;
import net.minecraft.world.level.levelgen.WorldgenRandom;
import net.minecraft.world.level.levelgen.XoroshiroRandomSource;
import net.minecraft.world.level.levelgen.placement.PlacedFeature;
import net.minecraft.world.level.levelgen.synth.NormalNoise;
import net.minecraft.world.level.levelgen.synth.PerlinNoise;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

public class Geode implements ModInitializer {
  private static final String MOD_ID = "geode";
  private static final Logger LOGGER = LoggerFactory.getLogger(MOD_ID);

  private static final long SEED = 0xDEAD_BEEFL;

  private static final AtomicInteger GEODE_COUNT = new AtomicInteger(0);
  private static final AtomicInteger BUDDING_COUNT = new AtomicInteger(0);

  private static WorldGenLevel fakeWorld(MinecraftServer server) {
    return (WorldGenLevel)
        Proxy.newProxyInstance(
            WorldGenLevel.class.getClassLoader(),
            new Class<?>[] {WorldGenLevel.class},
            (proxy, method, args) -> {
              String name = method.getName();
              return switch (name) {
                case "getSeed" -> SEED;
                case "getBlockState" -> Blocks.STONE.defaultBlockState();
                case "setBlock" -> {
                  BlockState state = (BlockState) args[1];
                  if (state.is(Blocks.BUDDING_AMETHYST)) {
                    BUDDING_COUNT.incrementAndGet();
                  }
                  yield true;
                }
                default -> method.invoke(server.overworld(), args);
              };
            });
  }

  private static void testRandom() {
    LOGGER.info("Random:");
    WorldgenRandom random = new WorldgenRandom(new LegacyRandomSource(SEED));
    LOGGER.info("  nextInt power-of-2: {}", random.nextInt(256));
    LOGGER.info("  nextInt non power-of-2: {}", random.nextInt(192));
    LOGGER.info("  nextBetween: {}", Mth.randomBetweenInclusive(random, 192, 256));
    LOGGER.info("  nextLong: {}", random.nextLong());
    LOGGER.info("  nextFloat: {}", (random.nextFloat()));
    LOGGER.info("  nextDouble: {}", (random.nextDouble()));
    random.consumeCount(256);
    LOGGER.info("  skip: {}", random.nextInt(256));
  }

  private static void testJavaRandom() {
    LOGGER.info("Java Random:");
    WorldgenRandom random = new WorldgenRandom(new LegacyRandomSource(SEED));

    int bits = 0;
    for (int i = 0; i < 256; i++) {
      bits = random.next(32);
    }
    LOGGER.info("  random progression: {}", bits);

    random.setSeed(SEED);
    RandomSource forkedRandom = random.forkPositional().fromHashOf("octave_-4");
    LOGGER.info("  fork: {}", forkedRandom.nextLong());
  }

  private static void testXoroshiro128Random() {
    LOGGER.info("Xoroshiro128++:");
    WorldgenRandom random = new WorldgenRandom(new XoroshiroRandomSource(SEED));

    int bits = 0;
    for (int i = 0; i < 256; i++) {
      bits = random.next(32);
    }
    LOGGER.info("  random progression: {}", bits);
  }

  private static void testPerlinNoise() {
    PerlinNoise perlin =
        PerlinNoise.createLegacyForLegacyNetherBiome(
            new WorldgenRandom(new LegacyRandomSource(SEED)), 0, DoubleList.of(1.0));
    LOGGER.info("Perlin Noise:");
    LOGGER.info("  0: {}", perlin.getValue(0.0, 0.0, 0.0));
    LOGGER.info("  -half: {}", perlin.getValue(-0.5, -0.5, -0.5));
    LOGGER.info("  wb: {}", perlin.getValue(-29_999_999.5, 64.25, 29_999_999.125));
  }

  private static void testNormalNoise() {
    int octave = -4;
    double amplitude = 1.0;

    // undo scaling done by noise
    double frequency = Math.pow(2.0, octave);
    double input_factor = 1.0181268882175227;
    double scale = 1.0 / (frequency * input_factor);

    NormalNoise noise =
        NormalNoise.create(new WorldgenRandom(new LegacyRandomSource(SEED)), octave, amplitude);
    LOGGER.info("Normal Noise:");
    LOGGER.info("  0: {}", noise.getValue(0.0, 0.0, 0.0));
    LOGGER.info("  -1: {}", noise.getValue(-scale, -scale, -scale));
    LOGGER.info("  mix: {}", noise.getValue(0.5 * scale, 0.25 * scale, 0.125 * scale));
  }

  private static <T> Registry<T> lookupRegistry(MinecraftServer server, String registryName) {
    ResourceKey<Registry<T>> registryKey =
        ResourceKey.createRegistryKey(new ResourceLocation(registryName));
    return server.registryAccess().registryOrThrow(registryKey);
  }

  private static <T> T lookup(MinecraftServer server, String registryName, String resourceName) {
    Registry<T> registry = lookupRegistry(server, registryName);
    return registry.get(new ResourceLocation(resourceName));
  }

  private static int[] getSalt(MinecraftServer server) {
    List<StepFeatureData> featuresByStep =
        server.getWorldData().worldGenSettings().overworld().getBiomeSource().featuresPerStep();
    ResourceLocation geodeKey = ResourceLocation.tryParse("amethyst_geode");

    int[] salt = new int[2];
    for (StepFeatureData genStep : featuresByStep) {
      salt[1] = 0;
      for (PlacedFeature placedFeature : genStep.features()) {
        if (placedFeature.feature().is(geodeKey)) {
          LOGGER.info("  salt: {}", 10000 * salt[0] + salt[1]);
          return salt;
        }
        salt[1]++;
      }
      salt[0]++;
    }

    return salt;
  }

  private static void testGeode(MinecraftServer server) {
    LOGGER.info("Geode:");
    WorldGenLevel world = fakeWorld(server);
    WorldgenRandom random = new WorldgenRandom(new XoroshiroRandomSource(SEED));
    ChunkGenerator generator = server.getWorldData().worldGenSettings().overworld();
    PlacedFeature geode = lookup(server, "worldgen/placed_feature", "amethyst_geode");

    int[] salt = getSalt(server);
    int range = 32;
    for (int z = -range; z <= range; z++) {
      for (int x = -range; x <= range; x++) {
        BlockPos blockPos = new BlockPos(x * 16, 0, z * 16);
        long l = random.setDecorationSeed(SEED, blockPos.getX(), blockPos.getZ());
        random.setFeatureSeed(l, salt[1], salt[0]);
        if (geode.placeWithBiomeCheck(world, generator, random, blockPos)) {
          GEODE_COUNT.incrementAndGet();
        }
      }
    }

    LOGGER.info("  geode count: {}", GEODE_COUNT);
    LOGGER.info("  budding count: {}", BUDDING_COUNT);
  }

  @Override
  public void onInitialize() {
    ServerLifecycleEvents.SERVER_STARTED.register(
        server -> {
          LOGGER.info("=== Testing for geode-finder:");
          testRandom();
          testJavaRandom();
          testXoroshiro128Random();
          testPerlinNoise();
          testNormalNoise();
          testGeode(server);
        });
  }
}
