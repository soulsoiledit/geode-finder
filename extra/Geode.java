// INFO: Configure as the entrypoint in a Fabric project to verify whether new MC versions will need to be added to geode-finder. May need modifications.

import net.fabricmc.api.ModInitializer;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import net.minecraft.util.RandomSource;
import net.minecraft.world.level.levelgen.LegacyRandomSource;
import net.minecraft.world.level.levelgen.WorldgenRandom;
import net.minecraft.world.level.levelgen.XoroshiroRandomSource;

import java.util.function.Function;
import net.minecraft.world.level.levelgen.synth.ImprovedNoise;
import net.minecraft.world.level.levelgen.synth.NormalNoise;
import net.minecraft.world.level.levelgen.synth.PerlinNoise;

import java.lang.reflect.Proxy;
import net.fabricmc.fabric.api.event.lifecycle.v1.ServerLifecycleEvents;
import net.minecraft.server.MinecraftServer;
import net.minecraft.world.level.WorldGenLevel;
import net.minecraft.world.level.block.Blocks;
import net.minecraft.world.level.block.state.BlockState;

import net.minecraft.core.Holder.Reference;
import net.minecraft.core.Registry;
import net.minecraft.core.RegistryAccess;
import net.minecraft.resources.Identifier;
import net.minecraft.resources.ResourceKey;

import java.util.List;
import net.minecraft.core.Holder;
import net.minecraft.core.HolderSet;
import net.minecraft.world.level.biome.Biome;
import net.minecraft.world.level.levelgen.placement.PlacedFeature;

import java.util.concurrent.atomic.AtomicInteger;
import net.minecraft.core.BlockPos;
import net.minecraft.world.level.chunk.ChunkGenerator;
import net.minecraft.world.level.dimension.LevelStem;
import net.minecraft.world.level.levelgen.NoiseGeneratorSettings;

public class Geode implements ModInitializer {
	private static final String MOD_ID = "geode";
	private static final Logger LOGGER = LoggerFactory.getLogger(MOD_ID);

	private static final long SEED = 0;
	private static final int NEAR = 0;
	private static final int FAR = 1875000;
	private static final int OFFSET = -4;
	private static final double AMPLITUDES = 1.0;
	private static final int RANGE = 16;

	private static final AtomicInteger GEODE_COUNT = new AtomicInteger(0);
	private static final AtomicInteger BUDDING_COUNT = new AtomicInteger(0);

	private static WorldGenLevel fakeWorld(MinecraftServer server) {
		return (WorldGenLevel) Proxy.newProxyInstance(WorldGenLevel.class.getClassLoader(),
				new Class<?>[] { WorldGenLevel.class }, (proxy, method, args) -> {
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

	private static void testRandom(String name, RandomSource r, boolean fork) {
		WorldgenRandom random = new WorldgenRandom(r);
		LOGGER.info("--- Random {}:", name);
		LOGGER.info("nextBits: {}", random.next(32));
		LOGGER.info("nextInt: {}", random.nextInt(256));
		LOGGER.info("nextBetween: {}", random.nextIntBetweenInclusive(16, 64));
		LOGGER.info("nextLong: {}", random.nextLong());
		LOGGER.info("nextFloat: {}", (random.nextFloat()));
		LOGGER.info("nextDouble: {}", (random.nextDouble()));
		random.consumeCount(1);
		LOGGER.info("skip nextInt: {}", random.nextInt(256));
		if (fork) {
			RandomSource forkedRandom = random.forkPositional().fromHashOf("octave_-4");
			LOGGER.info("forked nextInt: {}", forkedRandom.nextInt(256));
		}
	}

	@FunctionalInterface
	private static interface Noise {
		double getValue(double x, double y, double z);
	}

	private static void testNoise(String name, Function<WorldgenRandom, Noise> noiseProvider) {
		double zero = 0.0;
		double min = -30000000;
		double max = 30000000;

		Noise noise = noiseProvider.apply(new WorldgenRandom(new LegacyRandomSource(SEED)));
		LOGGER.info("--- Noise {}:", name);
		LOGGER.info("0: {}", noise.getValue(zero, zero, zero));
		LOGGER.info("min: {}", noise.getValue(min, min, min));
		LOGGER.info("max: {}", noise.getValue(max, max, max));
		LOGGER.info("all: {}", noise.getValue(min, zero, max));
	}

	private static <T> Registry<T> lookup(MinecraftServer server, String registryName) {
		RegistryAccess vanilla = server.registryAccess();
		ResourceKey<Registry<T>> registryKey = ResourceKey.createRegistryKey(Identifier.parse(registryName));
		return vanilla.lookupOrThrow(registryKey);
	}

	private static <T> T lookup(MinecraftServer server, String registryName, String resourceName) {
		RegistryAccess vanilla = server.registryAccess();
		ResourceKey<Registry<T>> registryKey = ResourceKey.createRegistryKey(Identifier.parse(registryName));
		Registry<T> registry = vanilla.lookupOrThrow(registryKey);

		Reference<T> holder = registry.get(Identifier.parse(resourceName)).orElseThrow();
		return holder.value();
	}

	private static int[] getSalt(MinecraftServer server) {
		Biome plains = lookup(server, "worldgen/biome", "plains");
		List<HolderSet<PlacedFeature>> featuresByStep = plains.getGenerationSettings().features();
		Identifier geodeLocation = Identifier.parse("amethyst_geode");

		int[] salt = new int[2];
		for (HolderSet<PlacedFeature> genStep : featuresByStep) {
			for (Holder<PlacedFeature> feature : genStep) {
				if (feature.is(geodeLocation)) {
					LOGGER.info("salt: {}", 10000 * salt[0] + salt[1]);
					return salt;
				}
				salt[1]++;
			}
			salt[0]++;
		}

		return salt;
	}

	private static void testGeode(String name, int center, MinecraftServer server, int[] salt) {
		LOGGER.info("--- Geode {}:", name);
		GEODE_COUNT.set(0);
		BUDDING_COUNT.set(0);

		WorldGenLevel world = fakeWorld(server);
		ChunkGenerator generator = ((LevelStem) lookup(server, "dimension", "overworld")).generator();
		NoiseGeneratorSettings owSettings = lookup(server, "worldgen/noise_settings", "overworld");
		WorldgenRandom random = new WorldgenRandom(owSettings.getRandomSource().newInstance(SEED));

		PlacedFeature geode = lookup(server, "worldgen/placed_feature", "amethyst_geode");
		Registry<PlacedFeature> placedRegistry = lookup(server, "worldgen/placed_feature");

		int min = center - RANGE;
		int max = center + RANGE;
		for (int i = min; i <= max; i++) {
			for (int j = min; j <= max; j++) {
				BlockPos blockPos = new BlockPos(i << 4, 10, j << 4);
				long l = random.setDecorationSeed(SEED, blockPos.getX(), blockPos.getZ());
				random.setFeatureSeed(l, salt[1], salt[0]);
				if (geode.placeWithBiomeCheck(world, generator, random, blockPos)) {
					GEODE_COUNT.incrementAndGet();
				}
			}
		}

		LOGGER.info("geode count: {}", GEODE_COUNT);
		LOGGER.info("amethyst count: {}", BUDDING_COUNT);
	}

	@Override
	public void onInitialize() {
		ServerLifecycleEvents.SERVER_STARTED.register(server -> {
			LOGGER.info("Testing for geode-finder:");
			testRandom("xoroshiro128++", new XoroshiroRandomSource(SEED), false);
			testRandom("java", new LegacyRandomSource(SEED), true);
			testNoise("improved", r -> new ImprovedNoise(r)::noise);
			testNoise("perlin", r -> PerlinNoise.create(r, OFFSET, AMPLITUDES)::getValue);
			testNoise("normal", r -> NormalNoise.create(r, OFFSET, AMPLITUDES)::getValue);
			int[] salt = getSalt(server);
			testGeode("near", NEAR, server, salt);
			testGeode("near", FAR, server, salt);
		});
	}
}
