
/*
    pbrt source code is Copyright(c) 1998-2016
                        Matt Pharr, Greg Humphreys, and Wenzel Jakob.

    This file is part of pbrt.

    Redistribution and use in source and binary forms, with or without
    modification, are permitted provided that the following conditions are
    met:

    - Redistributions of source code must retain the above copyright
      notice, this list of conditions and the following disclaimer.

    - Redistributions in binary form must reproduce the above copyright
      notice, this list of conditions and the following disclaimer in the
      documentation and/or other materials provided with the distribution.

    THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
    IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
    TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
    PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
    HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
    SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
    LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
    DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
    THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
    (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
    OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

 */

// core/api.cpp*
//
// 此文件是 pbrt 渲染器的核心 API 实现层。
// 实现了场景描述语言中所有 pbrt* 前缀的命令，包括：
// 变换操作（Translate/Rotate/Scale 等）、相机/胶片/采样器/积分器等
// 组件的创建与配置、材质系统、灯光系统以及场景构建与渲染控制。
// 同时也包含 ObjectInstance、AttributeBegin/End 等层次化场景管理功能。
//

#include "api.h"
#include "parallel.h"
#include "paramset.h"
#include "spectrum.h"
#include "scene.h"
#include "film.h"
#include "medium.h"
#include "stats.h"

// API Additional Headers
#include "accelerators/bvh.h"
#include "accelerators/kdtreeaccel.h"
#include "cameras/environment.h"
#include "cameras/orthographic.h"
#include "cameras/perspective.h"
#include "cameras/realistic.h"
#include "filters/box.h"
#include "filters/gaussian.h"
#include "filters/mitchell.h"
#include "filters/sinc.h"
#include "filters/triangle.h"
#include "integrators/bdpt.h"
#include "integrators/directlighting.h"
#include "integrators/mlt.h"
#include "integrators/ao.h"
#include "integrators/path.h"
#include "integrators/sppm.h"
#include "integrators/volpath.h"
#include "integrators/whitted.h"
#include "lights/diffuse.h"
#include "lights/distant.h"
#include "lights/goniometric.h"
#include "lights/infinite.h"
#include "lights/point.h"
#include "lights/projection.h"
#include "lights/spot.h"
#include "materials/disney.h"
#include "materials/fourier.h"
#include "materials/glass.h"
#include "materials/hair.h"
#include "materials/kdsubsurface.h"
#include "materials/matte.h"
#include "materials/metal.h"
#include "materials/mirror.h"
#include "materials/mixmat.h"
#include "materials/plastic.h"
#include "materials/substrate.h"
#include "materials/subsurface.h"
#include "materials/translucent.h"
#include "materials/uber.h"
#include "samplers/halton.h"
#include "samplers/maxmin.h"
#include "samplers/random.h"
#include "samplers/sobol.h"
#include "samplers/stratified.h"
#include "samplers/zerotwosequence.h"
#include "shapes/cone.h"
#include "shapes/curve.h"
#include "shapes/cylinder.h"
#include "shapes/disk.h"
#include "shapes/heightfield.h"
#include "shapes/hyperboloid.h"
#include "shapes/loopsubdiv.h"
#include "shapes/nurbs.h"
#include "shapes/paraboloid.h"
#include "shapes/sphere.h"
#include "shapes/triangle.h"
#include "shapes/plymesh.h"
#include "textures/bilerp.h"
#include "textures/checkerboard.h"
#include "textures/constant.h"
#include "textures/dots.h"
#include "textures/fbm.h"
#include "textures/imagemap.h"
#include "textures/marble.h"
#include "textures/mix.h"
#include "textures/ptex.h"
#include "textures/scale.h"
#include "textures/uv.h"
#include "textures/windy.h"
#include "textures/wrinkled.h"
#include "media/grid.h"
#include "media/homogeneous.h"

#include <map>
#include <stdio.h>

namespace pbrt {

// API Global Variables
Options PbrtOptions;

// API Local Classes
PBRT_CONSTEXPR int MaxTransforms = 2;
PBRT_CONSTEXPR int StartTransformBits = 1 << 0;
PBRT_CONSTEXPR int EndTransformBits = 1 << 1;
PBRT_CONSTEXPR int AllTransformsBits = (1 << MaxTransforms) - 1;
struct TransformSet {
    // TransformSet Public Methods
    Transform &operator[](int i) {
        CHECK_GE(i, 0);
        CHECK_LT(i, MaxTransforms);
        return t[i];
    }
    const Transform &operator[](int i) const {
        CHECK_GE(i, 0);
        CHECK_LT(i, MaxTransforms);
        return t[i];
    }
    friend TransformSet Inverse(const TransformSet &ts) {
        TransformSet tInv;
        for (int i = 0; i < MaxTransforms; ++i) tInv.t[i] = Inverse(ts.t[i]);
        return tInv;
    }
    bool IsAnimated() const {
        for (int i = 0; i < MaxTransforms - 1; ++i)
            if (t[i] != t[i + 1]) return true;
        return false;
    }

  private:
    Transform t[MaxTransforms];
};

struct RenderOptions {
    // RenderOptions Public Methods
    Integrator *MakeIntegrator() const;
    Scene *MakeScene();
    Camera *MakeCamera() const;

    // RenderOptions Public Data
    Float transformStartTime = 0, transformEndTime = 1;
    std::string FilterName = "box";
    ParamSet FilterParams;
    std::string FilmName = "image";
    ParamSet FilmParams;
    std::string SamplerName = "halton";
    ParamSet SamplerParams;
    std::string AcceleratorName = "bvh";
    ParamSet AcceleratorParams;
    std::string IntegratorName = "path";
    ParamSet IntegratorParams;
    std::string CameraName = "perspective";
    ParamSet CameraParams;
    TransformSet CameraToWorld;
    std::map<std::string, std::shared_ptr<Medium>> namedMedia;
    std::vector<std::shared_ptr<Light>> lights;
    std::vector<std::shared_ptr<Primitive>> primitives;
    std::map<std::string, std::vector<std::shared_ptr<Primitive>>> instances;
    std::vector<std::shared_ptr<Primitive>> *currentInstance = nullptr;
    bool haveScatteringMedia = false;
};

// MaterialInstance represents both an instance of a material as well as
// the information required to create another instance of it (possibly with
// different parameters from the shape).
struct MaterialInstance {
    MaterialInstance() = default;
    MaterialInstance(const std::string &name, const std::shared_ptr<Material> &mtl,
                     ParamSet params)
        : name(name), material(mtl), params(std::move(params)) {}

    std::string name;
    std::shared_ptr<Material> material;
    ParamSet params;
};

struct GraphicsState {
    // Graphics State Methods
    GraphicsState()
        : floatTextures(std::make_shared<FloatTextureMap>()),
          spectrumTextures(std::make_shared<SpectrumTextureMap>()),
          namedMaterials(std::make_shared<NamedMaterialMap>()) {
        ParamSet empty;
        TextureParams tp(empty, empty, *floatTextures, *spectrumTextures);
        std::shared_ptr<Material> mtl(CreateMatteMaterial(tp));
        currentMaterial = std::make_shared<MaterialInstance>("matte", mtl, ParamSet());
    }
    std::shared_ptr<Material> GetMaterialForShape(const ParamSet &geomParams);
    MediumInterface CreateMediumInterface();

    // Graphics State
    std::string currentInsideMedium, currentOutsideMedium;

    // Updated after book publication: floatTextures, spectrumTextures, and
    // namedMaterials are all implemented using a "copy on write" approach
    // for more efficient GraphicsState management.  When state is pushed
    // in pbrtAttributeBegin(), we don't immediately make a copy of these
    // maps, but instead record that each one is shared.  Only if an item
    // is added to one is a unique copy actually made.
    using FloatTextureMap = std::map<std::string, std::shared_ptr<Texture<Float>>>;
    std::shared_ptr<FloatTextureMap> floatTextures;
    bool floatTexturesShared = false;

    using SpectrumTextureMap = std::map<std::string, std::shared_ptr<Texture<Spectrum>>>;
    std::shared_ptr<SpectrumTextureMap> spectrumTextures;
    bool spectrumTexturesShared = false;

    using NamedMaterialMap = std::map<std::string, std::shared_ptr<MaterialInstance>>;
    std::shared_ptr<NamedMaterialMap> namedMaterials;
    bool namedMaterialsShared = false;

    std::shared_ptr<MaterialInstance> currentMaterial;
    ParamSet areaLightParams;
    std::string areaLight;
    bool reverseOrientation = false;
};

STAT_MEMORY_COUNTER("Memory/TransformCache", transformCacheBytes);
STAT_PERCENT("Scene/TransformCache hits", nTransformCacheHits, nTransformCacheLookups);
STAT_INT_DISTRIBUTION("Scene/Probes per TransformCache lookup", transformCacheProbes);

// Note: TransformCache has been reimplemented and has a slightly different
// interface compared to the version described in the third edition of
// Physically Based Rendering.  The new version is more efficient in both
// space and memory, which is helpful for highly complex scenes.
//
// The new implementation uses a hash table to store Transforms (rather
// than a std::map, which generally uses a red-black tree).  Further,
// it doesn't always store the inverse of the transform; if a caller
// wants the inverse as well, they are responsible for storing it.
//
// The hash table size is always a power of two, allowing for the use of a
// bitwise AND to turn hash values into table offsets.  Quadratic probing
// is used when there is a hash collision.
class TransformCache {
  public:
    TransformCache()
        : hashTable(512), hashTableOccupancy(0) {}

    // TransformCache Public Methods
    Transform *Lookup(const Transform &t) {
        ++nTransformCacheLookups;

        int offset = Hash(t) & (hashTable.size() - 1);
        int step = 1;
        while (true) {
            // Keep looking until we find the Transform or determine that
            // it's not present.
            if (!hashTable[offset] || *hashTable[offset] == t)
                break;
            // Advance using quadratic probing.
            offset = (offset + step * step) & (hashTable.size() - 1);
            ++step;
        }
        ReportValue(transformCacheProbes, step);
        Transform *tCached = hashTable[offset];
        if (tCached)
            ++nTransformCacheHits;
        else {
            tCached = arena.Alloc<Transform>();
            *tCached = t;
            Insert(tCached);
        }
        return tCached;
    }

    void Clear() {
        transformCacheBytes += arena.TotalAllocated() + hashTable.size() * sizeof(Transform *);
        hashTable.clear();
        hashTable.resize(512);
        hashTableOccupancy = 0;
        arena.Reset();
    }

  private:
    void Insert(Transform *tNew);
    void Grow();

    static uint64_t Hash(const Transform &t) {
        const char *ptr = (const char *)(&t.GetMatrix());
        size_t size = sizeof(Matrix4x4);
        uint64_t hash = 14695981039346656037ull;
        while (size > 0) {
            hash ^= *ptr;
            hash *= 1099511628211ull;
            ++ptr;
            --size;
        }
        return hash;
    }

    // TransformCache Private Data
    std::vector<Transform *> hashTable;
    int hashTableOccupancy;
    MemoryArena arena;
};

void TransformCache::Insert(Transform *tNew) {
    if (++hashTableOccupancy == hashTable.size() / 2)
        Grow();

    int baseOffset = Hash(*tNew) & (hashTable.size() - 1);
    for (int nProbes = 0;; ++nProbes) {
        // Quadratic probing.
        int offset = (baseOffset + nProbes/2 + nProbes*nProbes/2) & (hashTable.size() - 1);
        if (hashTable[offset] == nullptr) {
            hashTable[offset] = tNew;
            return;
        }
    }
}

void TransformCache::Grow() {
    std::vector<Transform *> newTable(2 * hashTable.size());
    LOG(INFO) << "Growing transform cache hash table to " << newTable.size();

    // Insert current elements into newTable.
    for (Transform *tEntry : hashTable) {
        if (!tEntry) continue;

        int baseOffset = Hash(*tEntry) & (hashTable.size() - 1);
        for (int nProbes = 0;; ++nProbes) {
            // Quadratic probing.
            int offset = (baseOffset + nProbes/2 + nProbes*nProbes/2) & (hashTable.size() - 1);
            if (newTable[offset] == nullptr) {
                newTable[offset] = tEntry;
                break;
            }
        }
    }

    std::swap(hashTable, newTable);
}


// API Static Data
enum class APIState { Uninitialized, OptionsBlock, WorldBlock };
static APIState currentApiState = APIState::Uninitialized;
static TransformSet curTransform;
static uint32_t activeTransformBits = AllTransformsBits;
static std::map<std::string, TransformSet> namedCoordinateSystems;
static std::unique_ptr<RenderOptions> renderOptions;
static GraphicsState graphicsState;
static std::vector<GraphicsState> pushedGraphicsStates;
static std::vector<TransformSet> pushedTransforms;
static std::vector<uint32_t> pushedActiveTransformBits;
static TransformCache transformCache;
int catIndentCount = 0;

// API Forward Declarations
std::vector<std::shared_ptr<Shape>> MakeShapes(const std::string &name,
                                               const Transform *ObjectToWorld,
                                               const Transform *WorldToObject,
                                               bool reverseOrientation,
                                               const ParamSet &paramSet);

// API Macros
#define VERIFY_INITIALIZED(func)                           \
    if (!(PbrtOptions.cat || PbrtOptions.toPly) &&           \
        currentApiState == APIState::Uninitialized) {        \
        Error(                                             \
            "pbrtInit() must be before calling \"%s()\". " \
            "Ignoring.",                                   \
            func);                                         \
        return;                                            \
    } else /* swallow trailing semicolon */
#define VERIFY_OPTIONS(func)                             \
    VERIFY_INITIALIZED(func);                            \
    if (!(PbrtOptions.cat || PbrtOptions.toPly) &&       \
        currentApiState == APIState::WorldBlock) {       \
        Error(                                           \
            "Options cannot be set inside world block; " \
            "\"%s\" not allowed.  Ignoring.",            \
            func);                                       \
        return;                                          \
    } else /* swallow trailing semicolon */
#define VERIFY_WORLD(func)                                   \
    VERIFY_INITIALIZED(func);                                \
    if (!(PbrtOptions.cat || PbrtOptions.toPly) &&           \
        currentApiState == APIState::OptionsBlock) {         \
        Error(                                               \
            "Scene description must be inside world block; " \
            "\"%s\" not allowed. Ignoring.",                 \
            func);                                           \
        return;                                              \
    } else /* swallow trailing semicolon */
#define FOR_ACTIVE_TRANSFORMS(expr)           \
    for (int i = 0; i < MaxTransforms; ++i)   \
        if (activeTransformBits & (1 << i)) { \
            expr                              \
        }
#define WARN_IF_ANIMATED_TRANSFORM(func)                             \
    do {                                                             \
        if (curTransform.IsAnimated())                               \
            Warning(                                                 \
                "Animated transformations set; ignoring for \"%s\" " \
                "and using the start transform only",                \
                func);                                               \
    } while (false) /* swallow trailing semicolon */

// Object Creation Function Definitions

// 根据名称和参数集创建形状对象的工厂函数。
// 支持 sphere、cylinder、disk、cone、paraboloid、hyperboloid、curve、
// trianglemesh、plymesh、heightfield、loopsubdiv、nurbs 等形状类型。
// 同时处理 -toply 模式下的 PLY 文件导出。
std::vector<std::shared_ptr<Shape>> MakeShapes(const std::string &name,
                                               const Transform *object2world,
                                               const Transform *world2object,
                                               bool reverseOrientation,
                                               const ParamSet &paramSet) {
    std::vector<std::shared_ptr<Shape>> shapes;
    std::shared_ptr<Shape> s;
    if (name == "sphere")
        s = CreateSphereShape(object2world, world2object, reverseOrientation,
                              paramSet);
    // Create remaining single _Shape_ types
    else if (name == "cylinder")
        s = CreateCylinderShape(object2world, world2object, reverseOrientation,
                                paramSet);
    else if (name == "disk")
        s = CreateDiskShape(object2world, world2object, reverseOrientation,
                            paramSet);
    else if (name == "cone")
        s = CreateConeShape(object2world, world2object, reverseOrientation,
                            paramSet);
    else if (name == "paraboloid")
        s = CreateParaboloidShape(object2world, world2object,
                                  reverseOrientation, paramSet);
    else if (name == "hyperboloid")
        s = CreateHyperboloidShape(object2world, world2object,
                                   reverseOrientation, paramSet);
    if (s != nullptr) shapes.push_back(s);

    // Create multiple-_Shape_ types
    else if (name == "curve")
        shapes = CreateCurveShape(object2world, world2object,
                                  reverseOrientation, paramSet);
    else if (name == "trianglemesh") {
        if (PbrtOptions.toPly) {
            int nvi;
            const int *vi = paramSet.FindInt("indices", &nvi);

            if (nvi < 500) {
                // It's a small mesh; don't bother with a PLY file after all.
                printf("%*sShape \"%s\" ", catIndentCount, "", name.c_str());
                paramSet.Print(catIndentCount);
                printf("\n");
            } else {
                static int count = 1;
                const char *plyPrefix =
                    getenv("PLY_PREFIX") ? getenv("PLY_PREFIX") : "mesh";
                std::string fn = StringPrintf("%s_%05d.ply", plyPrefix, count++);

                int npi, nuvi, nsi, nni;
                const Point3f *P = paramSet.FindPoint3f("P", &npi);
                const Point2f *uvs = paramSet.FindPoint2f("uv", &nuvi);
                if (!uvs) uvs = paramSet.FindPoint2f("st", &nuvi);
                std::vector<Point2f> tempUVs;
                if (!uvs) {
                    const Float *fuv = paramSet.FindFloat("uv", &nuvi);
                    if (!fuv) fuv = paramSet.FindFloat("st", &nuvi);
                    if (fuv) {
                        nuvi /= 2;
                        tempUVs.reserve(nuvi);
                        for (int i = 0; i < nuvi; ++i)
                            tempUVs.push_back(Point2f(fuv[2 * i], fuv[2 * i + 1]));
                        uvs = &tempUVs[0];
                    }
                }
                const Normal3f *N = paramSet.FindNormal3f("N", &nni);
                const Vector3f *S = paramSet.FindVector3f("S", &nsi);
                int nfi;
                const int *faceIndices = paramSet.FindInt("faceIndices", &nfi);
                if (faceIndices) CHECK_EQ(nfi, nvi / 3);

                if (!WritePlyFile(fn.c_str(), nvi / 3, vi, npi, P, S, N, uvs,
                                  faceIndices))
                    Error("Unable to write PLY file \"%s\"", fn.c_str());

                ParamSet ps = paramSet;
                ps.EraseInt("indices");
                ps.ErasePoint3f("P");
                ps.ErasePoint2f("uv");
                ps.ErasePoint2f("st");
                ps.EraseNormal3f("N");
                ps.EraseVector3f("S");
                ps.EraseInt("faceIndices");

                printf("%*sShape \"plymesh\" \"string filename\" \"%s\" ",
                       catIndentCount, "", fn.c_str());
                ps.Print(catIndentCount);
                printf("\n");
            }
        } else
            shapes = CreateTriangleMeshShape(object2world, world2object,
                                             reverseOrientation, paramSet,
                                             &*graphicsState.floatTextures);
    } else if (name == "plymesh")
        shapes = CreatePLYMesh(object2world, world2object, reverseOrientation,
                               paramSet, &*graphicsState.floatTextures);
    else if (name == "heightfield")
        shapes = CreateHeightfield(object2world, world2object,
                                   reverseOrientation, paramSet);
    else if (name == "loopsubdiv")
        shapes = CreateLoopSubdiv(object2world, world2object,
                                  reverseOrientation, paramSet);
    else if (name == "nurbs")
        shapes = CreateNURBS(object2world, world2object, reverseOrientation,
                             paramSet);
    else
        Warning("Shape \"%s\" unknown.", name.c_str());
    return shapes;
}

STAT_COUNTER("Scene/Materials created", nMaterialsCreated);

// 根据材质名称和纹理参数创建材质的工厂函数。
// 支持 matte、plastic、translucent、glass、mirror、hair、disney、
// mix、metal、substrate、uber、subsurface、kdsubsurface、fourier 等材质类型。
// 对于 mix 材质，会查找 namedMaterial 中定义的两个子材质。
std::shared_ptr<Material> MakeMaterial(const std::string &name,
                                       const TextureParams &mp) {
    Material *material = nullptr;
    if (name == "" || name == "none")
        return nullptr;
    else if (name == "matte")
        material = CreateMatteMaterial(mp);
    else if (name == "plastic")
        material = CreatePlasticMaterial(mp);
    else if (name == "translucent")
        material = CreateTranslucentMaterial(mp);
    else if (name == "glass")
        material = CreateGlassMaterial(mp);
    else if (name == "mirror")
        material = CreateMirrorMaterial(mp);
    else if (name == "hair")
        material = CreateHairMaterial(mp);
    else if (name == "disney")
        material = CreateDisneyMaterial(mp);
    else if (name == "mix") {
        std::string m1 = mp.FindString("namedmaterial1", "");
        std::string m2 = mp.FindString("namedmaterial2", "");
        std::shared_ptr<Material> mat1, mat2;
        if (graphicsState.namedMaterials->find(m1) ==
            graphicsState.namedMaterials->end()) {
            Error("Named material \"%s\" undefined.  Using \"matte\"",
                  m1.c_str());
            mat1 = MakeMaterial("matte", mp);
        } else
            mat1 = (*graphicsState.namedMaterials)[m1]->material;

        if (graphicsState.namedMaterials->find(m2) ==
            graphicsState.namedMaterials->end()) {
            Error("Named material \"%s\" undefined.  Using \"matte\"",
                  m2.c_str());
            mat2 = MakeMaterial("matte", mp);
        } else
            mat2 = (*graphicsState.namedMaterials)[m2]->material;

        material = CreateMixMaterial(mp, mat1, mat2);
    } else if (name == "metal")
        material = CreateMetalMaterial(mp);
    else if (name == "substrate")
        material = CreateSubstrateMaterial(mp);
    else if (name == "uber")
        material = CreateUberMaterial(mp);
    else if (name == "subsurface")
        material = CreateSubsurfaceMaterial(mp);
    else if (name == "kdsubsurface")
        material = CreateKdSubsurfaceMaterial(mp);
    else if (name == "fourier")
        material = CreateFourierMaterial(mp);
    else {
        Warning("Material \"%s\" unknown. Using \"matte\".", name.c_str());
        material = CreateMatteMaterial(mp);
    }

    if ((name == "subsurface" || name == "kdsubsurface") &&
        (renderOptions->IntegratorName != "path" &&
         (renderOptions->IntegratorName != "volpath")))
        Warning(
            "Subsurface scattering material \"%s\" used, but \"%s\" "
            "integrator doesn't support subsurface scattering. "
            "Use \"path\" or \"volpath\".",
            name.c_str(), renderOptions->IntegratorName.c_str());

    mp.ReportUnused();
    if (!material) Error("Unable to create material \"%s\"", name.c_str());
    else ++nMaterialsCreated;
    return std::shared_ptr<Material>(material);
}

// 创建 Float 纹理的工厂函数。支持 constant、scale、mix、bilerp、
// imagemap、uv、checkerboard、dots、fbm、wrinkled、marble、windy、ptex 等纹理类型。
std::shared_ptr<Texture<Float>> MakeFloatTexture(const std::string &name,
                                                 const Transform &tex2world,
                                                 const TextureParams &tp) {
    Texture<Float> *tex = nullptr;
    if (name == "constant")
        tex = CreateConstantFloatTexture(tex2world, tp);
    else if (name == "scale")
        tex = CreateScaleFloatTexture(tex2world, tp);
    else if (name == "mix")
        tex = CreateMixFloatTexture(tex2world, tp);
    else if (name == "bilerp")
        tex = CreateBilerpFloatTexture(tex2world, tp);
    else if (name == "imagemap")
        tex = CreateImageFloatTexture(tex2world, tp);
    else if (name == "uv")
        tex = CreateUVFloatTexture(tex2world, tp);
    else if (name == "checkerboard")
        tex = CreateCheckerboardFloatTexture(tex2world, tp);
    else if (name == "dots")
        tex = CreateDotsFloatTexture(tex2world, tp);
    else if (name == "fbm")
        tex = CreateFBmFloatTexture(tex2world, tp);
    else if (name == "wrinkled")
        tex = CreateWrinkledFloatTexture(tex2world, tp);
    else if (name == "marble")
        tex = CreateMarbleFloatTexture(tex2world, tp);
    else if (name == "windy")
        tex = CreateWindyFloatTexture(tex2world, tp);
    else if (name == "ptex")
        tex = CreatePtexFloatTexture(tex2world, tp);
    else
        Warning("Float texture \"%s\" unknown.", name.c_str());
    tp.ReportUnused();
    return std::shared_ptr<Texture<Float>>(tex);
}

// 创建 Spectrum 纹理的工厂函数。支持的纹理类型与 Float 纹理一一对应。
std::shared_ptr<Texture<Spectrum>> MakeSpectrumTexture(
    const std::string &name, const Transform &tex2world,
    const TextureParams &tp) {
    Texture<Spectrum> *tex = nullptr;
    if (name == "constant")
        tex = CreateConstantSpectrumTexture(tex2world, tp);
    else if (name == "scale")
        tex = CreateScaleSpectrumTexture(tex2world, tp);
    else if (name == "mix")
        tex = CreateMixSpectrumTexture(tex2world, tp);
    else if (name == "bilerp")
        tex = CreateBilerpSpectrumTexture(tex2world, tp);
    else if (name == "imagemap")
        tex = CreateImageSpectrumTexture(tex2world, tp);
    else if (name == "uv")
        tex = CreateUVSpectrumTexture(tex2world, tp);
    else if (name == "checkerboard")
        tex = CreateCheckerboardSpectrumTexture(tex2world, tp);
    else if (name == "dots")
        tex = CreateDotsSpectrumTexture(tex2world, tp);
    else if (name == "fbm")
        tex = CreateFBmSpectrumTexture(tex2world, tp);
    else if (name == "wrinkled")
        tex = CreateWrinkledSpectrumTexture(tex2world, tp);
    else if (name == "marble")
        tex = CreateMarbleSpectrumTexture(tex2world, tp);
    else if (name == "windy")
        tex = CreateWindySpectrumTexture(tex2world, tp);
    else if (name == "ptex")
        tex = CreatePtexSpectrumTexture(tex2world, tp);
    else
        Warning("Spectrum texture \"%s\" unknown.", name.c_str());
    tp.ReportUnused();
    return std::shared_ptr<Texture<Spectrum>>(tex);
}

// 创建介质的工厂函数。支持 homogeneous（均匀）和 heterogeneous（非均匀/网格密度）介质。
// 从参数集中读取散射系数 sigma_a 和 sigma_s、各向异性参数 g 以及缩放因子。
std::shared_ptr<Medium> MakeMedium(const std::string &name,
                                   const ParamSet &paramSet,
                                   const Transform &medium2world) {
    Float sig_a_rgb[3] = {.0011f, .0024f, .014f},
          sig_s_rgb[3] = {2.55f, 3.21f, 3.77f};
    Spectrum sig_a = Spectrum::FromRGB(sig_a_rgb),
             sig_s = Spectrum::FromRGB(sig_s_rgb);
    std::string preset = paramSet.FindOneString("preset", "");
    bool found = GetMediumScatteringProperties(preset, &sig_a, &sig_s);
    if (preset != "" && !found)
        Warning("Material preset \"%s\" not found.  Using defaults.",
                preset.c_str());
    Float scale = paramSet.FindOneFloat("scale", 1.f);
    Float g = paramSet.FindOneFloat("g", 0.0f);
    sig_a = paramSet.FindOneSpectrum("sigma_a", sig_a) * scale;
    sig_s = paramSet.FindOneSpectrum("sigma_s", sig_s) * scale;
    Medium *m = NULL;
    if (name == "homogeneous") {
        m = new HomogeneousMedium(sig_a, sig_s, g);
    } else if (name == "heterogeneous") {
        int nitems;
        const Float *data = paramSet.FindFloat("density", &nitems);
        if (!data) {
            Error("No \"density\" values provided for heterogeneous medium?");
            return NULL;
        }
        int nx = paramSet.FindOneInt("nx", 1);
        int ny = paramSet.FindOneInt("ny", 1);
        int nz = paramSet.FindOneInt("nz", 1);
        Point3f p0 = paramSet.FindOnePoint3f("p0", Point3f(0.f, 0.f, 0.f));
        Point3f p1 = paramSet.FindOnePoint3f("p1", Point3f(1.f, 1.f, 1.f));
        if (nitems != nx * ny * nz) {
            Error(
                "GridDensityMedium has %d density values; expected nx*ny*nz = "
                "%d",
                nitems, nx * ny * nz);
            return NULL;
        }
        Transform data2Medium = Translate(Vector3f(p0)) *
                                Scale(p1.x - p0.x, p1.y - p0.y, p1.z - p0.z);
        m = new GridDensityMedium(sig_a, sig_s, g, nx, ny, nz,
                                  medium2world * data2Medium, data);
    } else
        Warning("Medium \"%s\" unknown.", name.c_str());
    paramSet.ReportUnused();
    return std::shared_ptr<Medium>(m);
}

// 创建非面光源的工厂函数。支持 point、spot、goniometric、projection、distant、infinite 等类型。
std::shared_ptr<Light> MakeLight(const std::string &name,
                                 const ParamSet &paramSet,
                                 const Transform &light2world,
                                 const MediumInterface &mediumInterface) {
    std::shared_ptr<Light> light;
    if (name == "point")
        light =
            CreatePointLight(light2world, mediumInterface.outside, paramSet);
    else if (name == "spot")
        light = CreateSpotLight(light2world, mediumInterface.outside, paramSet);
    else if (name == "goniometric")
        light = CreateGoniometricLight(light2world, mediumInterface.outside,
                                       paramSet);
    else if (name == "projection")
        light = CreateProjectionLight(light2world, mediumInterface.outside,
                                      paramSet);
    else if (name == "distant")
        light = CreateDistantLight(light2world, paramSet);
    else if (name == "infinite" || name == "exinfinite")
        light = CreateInfiniteLight(light2world, paramSet);
    else
        Warning("Light \"%s\" unknown.", name.c_str());
    paramSet.ReportUnused();
    return light;
}

// 创建面光源（AreaLight）的工厂函数。面光源与几何形状关联，为形状赋予自发光属性。
std::shared_ptr<AreaLight> MakeAreaLight(const std::string &name,
                                         const Transform &light2world,
                                         const MediumInterface &mediumInterface,
                                         const ParamSet &paramSet,
                                         const std::shared_ptr<Shape> &shape) {
    std::shared_ptr<AreaLight> area;
    if (name == "area" || name == "diffuse")
        area = CreateDiffuseAreaLight(light2world, mediumInterface.outside,
                                      paramSet, shape);
    else
        Warning("Area light \"%s\" unknown.", name.c_str());
    paramSet.ReportUnused();
    return area;
}

// 创建加速结构的工厂函数。支持 bvh（层次包围盒）和 kdtree（Kd 树）加速器。
std::shared_ptr<Primitive> MakeAccelerator(
    const std::string &name,
    std::vector<std::shared_ptr<Primitive>> prims,
    const ParamSet &paramSet) {
    std::shared_ptr<Primitive> accel;
    if (name == "bvh")
        accel = CreateBVHAccelerator(std::move(prims), paramSet);
    else if (name == "kdtree")
        accel = CreateKdTreeAccelerator(std::move(prims), paramSet);
    else
        Warning("Accelerator \"%s\" unknown.", name.c_str());
    paramSet.ReportUnused();
    return accel;
}

// 创建相机对象的工厂函数。支持 perspective（透视）、orthographic（正交）、
// realistic（真实感）和 environment（环境）相机类型。
Camera *MakeCamera(const std::string &name, const ParamSet &paramSet,
                   const TransformSet &cam2worldSet, Float transformStart,
                   Float transformEnd, Film *film) {
    Camera *camera = nullptr;
    MediumInterface mediumInterface = graphicsState.CreateMediumInterface();
    static_assert(MaxTransforms == 2,
                  "TransformCache assumes only two transforms");
    Transform *cam2world[2] = {
        transformCache.Lookup(cam2worldSet[0]),
        transformCache.Lookup(cam2worldSet[1])
    };
    AnimatedTransform animatedCam2World(cam2world[0], transformStart,
                                        cam2world[1], transformEnd);
    if (name == "perspective")
        camera = CreatePerspectiveCamera(paramSet, animatedCam2World, film,
                                         mediumInterface.outside);
    else if (name == "orthographic")
        camera = CreateOrthographicCamera(paramSet, animatedCam2World, film,
                                          mediumInterface.outside);
    else if (name == "realistic")
        camera = CreateRealisticCamera(paramSet, animatedCam2World, film,
                                       mediumInterface.outside);
    else if (name == "environment")
        camera = CreateEnvironmentCamera(paramSet, animatedCam2World, film,
                                         mediumInterface.outside);
    else
        Warning("Camera \"%s\" unknown.", name.c_str());
    paramSet.ReportUnused();
    return camera;
}

// 创建采样器的工厂函数。支持 lowdiscrepancy、maxmindist、halton、
// sobol、random、stratified 等采样策略。
std::shared_ptr<Sampler> MakeSampler(const std::string &name,
                                     const ParamSet &paramSet,
                                     const Film *film) {
    Sampler *sampler = nullptr;
    if (name == "lowdiscrepancy" || name == "02sequence")
        sampler = CreateZeroTwoSequenceSampler(paramSet);
    else if (name == "maxmindist")
        sampler = CreateMaxMinDistSampler(paramSet);
    else if (name == "halton")
        sampler = CreateHaltonSampler(paramSet, film->GetSampleBounds());
    else if (name == "sobol")
        sampler = CreateSobolSampler(paramSet, film->GetSampleBounds());
    else if (name == "random")
        sampler = CreateRandomSampler(paramSet);
    else if (name == "stratified")
        sampler = CreateStratifiedSampler(paramSet);
    else
        Warning("Sampler \"%s\" unknown.", name.c_str());
    paramSet.ReportUnused();
    return std::shared_ptr<Sampler>(sampler);
}

// 创建图像重构滤波器的工厂函数。支持 box、gaussian、mitchell、sinc、triangle 滤波器。
std::unique_ptr<Filter> MakeFilter(const std::string &name,
                                   const ParamSet &paramSet) {
    Filter *filter = nullptr;
    if (name == "box")
        filter = CreateBoxFilter(paramSet);
    else if (name == "gaussian")
        filter = CreateGaussianFilter(paramSet);
    else if (name == "mitchell")
        filter = CreateMitchellFilter(paramSet);
    else if (name == "sinc")
        filter = CreateSincFilter(paramSet);
    else if (name == "triangle")
        filter = CreateTriangleFilter(paramSet);
    else {
        Error("Filter \"%s\" unknown.", name.c_str());
        exit(1);
    }
    paramSet.ReportUnused();
    return std::unique_ptr<Filter>(filter);
}

// 创建胶片对象的工厂函数。目前支持 "image" 类型胶片。
Film *MakeFilm(const std::string &name, const ParamSet &paramSet,
               std::unique_ptr<Filter> filter) {
    Film *film = nullptr;
    if (name == "image")
        film = CreateFilm(paramSet, std::move(filter));
    else
        Warning("Film \"%s\" unknown.", name.c_str());
    paramSet.ReportUnused();
    return film;
}

// API Function Definitions

// 初始化 pbrt 渲染系统。必须在所有 pbrt* API 调用之前调用。
// 设置渲染选项、初始化全局状态、采样光谱、并行系统和性能分析器。
void pbrtInit(const Options &opt) {
    PbrtOptions = opt;
    // API Initialization
    if (currentApiState != APIState::Uninitialized)
        Error("pbrtInit() has already been called.");
    currentApiState = APIState::OptionsBlock;
    renderOptions.reset(new RenderOptions);
    graphicsState = GraphicsState();
    catIndentCount = 0;

    // General \pbrt Initialization
    SampledSpectrum::Init();
    ParallelInit();  // Threads must be launched before the profiler is
                     // initialized.
    InitProfiler();
}

// 清理 pbrt 渲染系统。在渲染完成后调用，释放并行系统和性能分析器资源。
void pbrtCleanup() {
    // API Cleanup
    if (currentApiState == APIState::Uninitialized)
        Error("pbrtCleanup() called without pbrtInit().");
    else if (currentApiState == APIState::WorldBlock)
        Error("pbrtCleanup() called while inside world block.");
    currentApiState = APIState::Uninitialized;
    ParallelCleanup();
    CleanupProfiler();
}

// 将当前变换重置为单位矩阵。
void pbrtIdentity() {
    VERIFY_INITIALIZED("Identity");
    FOR_ACTIVE_TRANSFORMS(curTransform[i] = Transform();)
    if (PbrtOptions.cat || PbrtOptions.toPly)
        printf("%*sIdentity\n", catIndentCount, "");
}

// 在当前变换上附加一个平移变换。平移向量为 (dx, dy, dz)。
void pbrtTranslate(Float dx, Float dy, Float dz) {
    VERIFY_INITIALIZED("Translate");
    FOR_ACTIVE_TRANSFORMS(curTransform[i] = curTransform[i] *
                                            Translate(Vector3f(dx, dy, dz));)
    if (PbrtOptions.cat || PbrtOptions.toPly)
        printf("%*sTranslate %.9g %.9g %.9g\n", catIndentCount, "", dx, dy,
               dz);
}

// 用给定的 4x4 矩阵直接设置当前变换。
void pbrtTransform(Float tr[16]) {
    VERIFY_INITIALIZED("Transform");
    FOR_ACTIVE_TRANSFORMS(
        curTransform[i] = Transform(Matrix4x4(
            tr[0], tr[4], tr[8], tr[12], tr[1], tr[5], tr[9], tr[13], tr[2],
            tr[6], tr[10], tr[14], tr[3], tr[7], tr[11], tr[15]));)
    if (PbrtOptions.cat || PbrtOptions.toPly) {
        printf("%*sTransform [ ", catIndentCount, "");
        for (int i = 0; i < 16; ++i) printf("%.9g ", tr[i]);
        printf(" ]\n");
    }
}

// 在当前变换上右乘一个 4x4 矩阵。
void pbrtConcatTransform(Float tr[16]) {
    VERIFY_INITIALIZED("ConcatTransform");
    FOR_ACTIVE_TRANSFORMS(
        curTransform[i] =
            curTransform[i] *
            Transform(Matrix4x4(tr[0], tr[4], tr[8], tr[12], tr[1], tr[5],
                                tr[9], tr[13], tr[2], tr[6], tr[10], tr[14],
                                tr[3], tr[7], tr[11], tr[15]));)
    if (PbrtOptions.cat || PbrtOptions.toPly) {
        printf("%*sConcatTransform [ ", catIndentCount, "");
        for (int i = 0; i < 16; ++i) printf("%.9g ", tr[i]);
        printf(" ]\n");
    }
}

// 在当前变换上附加一个旋转变换。绕轴 (dx, dy, dz) 旋转 angle 度。
void pbrtRotate(Float angle, Float dx, Float dy, Float dz) {
    VERIFY_INITIALIZED("Rotate");
    FOR_ACTIVE_TRANSFORMS(curTransform[i] =
                              curTransform[i] *
                              Rotate(angle, Vector3f(dx, dy, dz));)
    if (PbrtOptions.cat || PbrtOptions.toPly)
        printf("%*sRotate %.9g %.9g %.9g %.9g\n", catIndentCount, "", angle,
               dx, dy, dz);
}

// 在当前变换上附加一个缩放变换。三个轴上的缩放因子分别为 sx, sy, sz。
void pbrtScale(Float sx, Float sy, Float sz) {
    VERIFY_INITIALIZED("Scale");
    FOR_ACTIVE_TRANSFORMS(curTransform[i] =
                              curTransform[i] * Scale(sx, sy, sz);)
    if (PbrtOptions.cat || PbrtOptions.toPly)
        printf("%*sScale %.9g %.9g %.9g\n", catIndentCount, "", sx, sy, sz);
}

// 设置观察变换。从视点 (ex,ey,ez) 看向目标点 (lx,ly,lz)，朝上方向为 (ux,uy,uz)。
void pbrtLookAt(Float ex, Float ey, Float ez, Float lx, Float ly, Float lz,
                Float ux, Float uy, Float uz) {
    VERIFY_INITIALIZED("LookAt");
    Transform lookAt =
        LookAt(Point3f(ex, ey, ez), Point3f(lx, ly, lz), Vector3f(ux, uy, uz));
    FOR_ACTIVE_TRANSFORMS(curTransform[i] = curTransform[i] * lookAt;);
    if (PbrtOptions.cat || PbrtOptions.toPly)
        printf(
            "%*sLookAt %.9g %.9g %.9g\n%*s%.9g %.9g %.9g\n"
            "%*s%.9g %.9g %.9g\n",
            catIndentCount, "", ex, ey, ez, catIndentCount + 8, "", lx, ly, lz,
            catIndentCount + 8, "", ux, uy, uz);
}

// 将当前变换保存为一个命名的坐标系，供后续 CoordSysTransform 引用。
void pbrtCoordinateSystem(const std::string &name) {
    VERIFY_INITIALIZED("CoordinateSystem");
    namedCoordinateSystems[name] = curTransform;
    if (PbrtOptions.cat || PbrtOptions.toPly)
        printf("%*sCoordinateSystem \"%s\"\n", catIndentCount, "",
               name.c_str());
}

// 将当前变换切换为之前通过 CoordinateSystem 保存的命名坐标系。
void pbrtCoordSysTransform(const std::string &name) {
    VERIFY_INITIALIZED("CoordSysTransform");
    if (namedCoordinateSystems.find(name) != namedCoordinateSystems.end())
        curTransform = namedCoordinateSystems[name];
    else
        Warning("Couldn't find named coordinate system \"%s\"", name.c_str());
    if (PbrtOptions.cat || PbrtOptions.toPly)
        printf("%*sCoordSysTransform \"%s\"\n", catIndentCount, "",
               name.c_str());
}

// 设置所有变换（开始时间和结束时间变换）均为活跃状态。
void pbrtActiveTransformAll() {
    activeTransformBits = AllTransformsBits;
    if (PbrtOptions.cat || PbrtOptions.toPly)
        printf("%*sActiveTransform All\n", catIndentCount, "");
}

// 设置只有结束时间变换为活跃状态。
void pbrtActiveTransformEndTime() {
    activeTransformBits = EndTransformBits;
    if (PbrtOptions.cat || PbrtOptions.toPly)
        printf("%*sActiveTransform EndTime\n", catIndentCount, "");
}

// 设置只有开始时间变换为活跃状态。
void pbrtActiveTransformStartTime() {
    activeTransformBits = StartTransformBits;
    if (PbrtOptions.cat || PbrtOptions.toPly)
        printf("%*sActiveTransform StartTime\n", catIndentCount, "");
}

// 设置动画变换的开始时间和结束时间。
void pbrtTransformTimes(Float start, Float end) {
    VERIFY_OPTIONS("TransformTimes");
    renderOptions->transformStartTime = start;
    renderOptions->transformEndTime = end;
    if (PbrtOptions.cat || PbrtOptions.toPly)
        printf("%*sTransformTimes %.9g %.9g\n", catIndentCount, "", start,
               end);
}

// 设置像素重构滤波器类型和参数（在 Options 块中调用）。
void pbrtPixelFilter(const std::string &name, const ParamSet &params) {
    VERIFY_OPTIONS("PixelFilter");
    renderOptions->FilterName = name;
    renderOptions->FilterParams = params;
    if (PbrtOptions.cat || PbrtOptions.toPly) {
        printf("%*sPixelFilter \"%s\" ", catIndentCount, "", name.c_str());
        params.Print(catIndentCount);
        printf("\n");
    }
}

// 设置胶片类型和参数（在 Options 块中调用）。
void pbrtFilm(const std::string &type, const ParamSet &params) {
    VERIFY_OPTIONS("Film");
    renderOptions->FilmParams = params;
    renderOptions->FilmName = type;
    if (PbrtOptions.cat || PbrtOptions.toPly) {
        printf("%*sFilm \"%s\" ", catIndentCount, "", type.c_str());
        params.Print(catIndentCount);
        printf("\n");
    }
}

// 设置采样器类型和参数（在 Options 块中调用）。
void pbrtSampler(const std::string &name, const ParamSet &params) {
    VERIFY_OPTIONS("Sampler");
    renderOptions->SamplerName = name;
    renderOptions->SamplerParams = params;
    if (PbrtOptions.cat || PbrtOptions.toPly) {
        printf("%*sSampler \"%s\" ", catIndentCount, "", name.c_str());
        params.Print(catIndentCount);
        printf("\n");
    }
}

// 设置加速结构类型和参数（在 Options 块中调用）。
void pbrtAccelerator(const std::string &name, const ParamSet &params) {
    VERIFY_OPTIONS("Accelerator");
    renderOptions->AcceleratorName = name;
    renderOptions->AcceleratorParams = params;
    if (PbrtOptions.cat || PbrtOptions.toPly) {
        printf("%*sAccelerator \"%s\" ", catIndentCount, "", name.c_str());
        params.Print(catIndentCount);
        printf("\n");
    }
}

// 设置积分器类型和参数（在 Options 块中调用）。
void pbrtIntegrator(const std::string &name, const ParamSet &params) {
    VERIFY_OPTIONS("Integrator");
    renderOptions->IntegratorName = name;
    renderOptions->IntegratorParams = params;
    if (PbrtOptions.cat || PbrtOptions.toPly) {
        printf("%*sIntegrator \"%s\" ", catIndentCount, "", name.c_str());
        params.Print(catIndentCount);
        printf("\n");
    }
}

// 设置相机类型和参数（在 Options 块中调用）。
// 同时将当前变换的逆矩阵设置为相机到世界的变换，并保存 "camera" 坐标系。
void pbrtCamera(const std::string &name, const ParamSet &params) {
    VERIFY_OPTIONS("Camera");
    renderOptions->CameraName = name;
    renderOptions->CameraParams = params;
    renderOptions->CameraToWorld = Inverse(curTransform);
    namedCoordinateSystems["camera"] = renderOptions->CameraToWorld;
    if (PbrtOptions.cat || PbrtOptions.toPly) {
        printf("%*sCamera \"%s\" ", catIndentCount, "", name.c_str());
        params.Print(catIndentCount);
        printf("\n");
    }
}

// 创建命名的介质对象，存储到 renderOptions 中供后续 MediumInterface 引用。
void pbrtMakeNamedMedium(const std::string &name, const ParamSet &params) {
    VERIFY_INITIALIZED("MakeNamedMedium");
    WARN_IF_ANIMATED_TRANSFORM("MakeNamedMedium");
    std::string type = params.FindOneString("type", "");
    if (type == "")
        Error("No parameter string \"type\" found in MakeNamedMedium");
    else {
        std::shared_ptr<Medium> medium =
            MakeMedium(type, params, curTransform[0]);
        if (medium) renderOptions->namedMedia[name] = medium;
    }
    if (PbrtOptions.cat || PbrtOptions.toPly) {
        printf("%*sMakeNamedMedium \"%s\" ", catIndentCount, "", name.c_str());
        params.Print(catIndentCount);
        printf("\n");
    }
}

// 设置当前图形状态中表面内部和外部介质的名称。
// 这些名称需与 MakeNamedMedium 中定义的名称匹配。
void pbrtMediumInterface(const std::string &insideName,
                         const std::string &outsideName) {
    VERIFY_INITIALIZED("MediumInterface");
    graphicsState.currentInsideMedium = insideName;
    graphicsState.currentOutsideMedium = outsideName;
    renderOptions->haveScatteringMedia = true;
    if (PbrtOptions.cat || PbrtOptions.toPly)
        printf("%*sMediumInterface \"%s\" \"%s\"\n", catIndentCount, "",
               insideName.c_str(), outsideName.c_str());
}

// 开始世界块（World Block）。在 WorldBegin 之后所有的变换和几何操作
// 都是在世界空间中进行的，直到 WorldEnd。
void pbrtWorldBegin() {
    VERIFY_OPTIONS("WorldBegin");
    currentApiState = APIState::WorldBlock;
    for (int i = 0; i < MaxTransforms; ++i) curTransform[i] = Transform();
    activeTransformBits = AllTransformsBits;
    namedCoordinateSystems["world"] = curTransform;
    if (PbrtOptions.cat || PbrtOptions.toPly)
        printf("\n\nWorldBegin\n\n");
}

// 开始属性块（Attribute Block）。保存当前图形状态、变换矩阵和活跃变换位，
// 推入栈中以便后续 AttributeEnd 恢复。
void pbrtAttributeBegin() {
    VERIFY_WORLD("AttributeBegin");
    pushedGraphicsStates.push_back(graphicsState);
    graphicsState.floatTexturesShared = graphicsState.spectrumTexturesShared =
        graphicsState.namedMaterialsShared = true;
    pushedTransforms.push_back(curTransform);
    pushedActiveTransformBits.push_back(activeTransformBits);
    if (PbrtOptions.cat || PbrtOptions.toPly) {
        printf("\n%*sAttributeBegin\n", catIndentCount, "");
        catIndentCount += 4;
    }
}

// 结束属性块。从栈中恢复之前保存的图形状态、变换矩阵和活跃变换位。
void pbrtAttributeEnd() {
    VERIFY_WORLD("AttributeEnd");
    if (!pushedGraphicsStates.size()) {
        Error(
            "Unmatched pbrtAttributeEnd() encountered. "
            "Ignoring it.");
        return;
    }
    graphicsState = std::move(pushedGraphicsStates.back());
    pushedGraphicsStates.pop_back();
    curTransform = pushedTransforms.back();
    pushedTransforms.pop_back();
    activeTransformBits = pushedActiveTransformBits.back();
    pushedActiveTransformBits.pop_back();
    if (PbrtOptions.cat || PbrtOptions.toPly) {
        catIndentCount -= 4;
        printf("%*sAttributeEnd\n", catIndentCount, "");
    }
}

// 开始变换块（Transform Block）。保存当前变换到栈中，
// 允许临时修改变换并在 TransformEnd 时恢复。
void pbrtTransformBegin() {
    VERIFY_WORLD("TransformBegin");
    pushedTransforms.push_back(curTransform);
    pushedActiveTransformBits.push_back(activeTransformBits);
    if (PbrtOptions.cat || PbrtOptions.toPly) {
        printf("%*sTransformBegin\n", catIndentCount, "");
        catIndentCount += 4;
    }
}

// 结束变换块。从栈中恢复之前保存的变换矩阵和活跃变换位。
void pbrtTransformEnd() {
    VERIFY_WORLD("TransformEnd");
    if (!pushedTransforms.size()) {
        Error(
            "Unmatched pbrtTransformEnd() encountered. "
            "Ignoring it.");
        return;
    }
    curTransform = pushedTransforms.back();
    pushedTransforms.pop_back();
    activeTransformBits = pushedActiveTransformBits.back();
    pushedActiveTransformBits.pop_back();
    if (PbrtOptions.cat || PbrtOptions.toPly) {
        catIndentCount -= 4;
        printf("%*sTransformEnd\n", catIndentCount, "");
    }
}

// 定义纹理对象。支持 float 和 color/spectrum 两种纹理类型。
// 纹理通过 MakeFloatTexture 或 MakeSpectrumTexture 创建，并存入当前图形状态的纹理映射表。
void pbrtTexture(const std::string &name, const std::string &type,
                 const std::string &texname, const ParamSet &params) {
    VERIFY_WORLD("Texture");
    if (PbrtOptions.cat || PbrtOptions.toPly) {
        printf("%*sTexture \"%s\" \"%s\" \"%s\" ", catIndentCount, "",
               name.c_str(), type.c_str(), texname.c_str());
        params.Print(catIndentCount);
        printf("\n");
        return;
    }

    TextureParams tp(params, params, *graphicsState.floatTextures,
                     *graphicsState.spectrumTextures);
    if (type == "float") {
        // Create _Float_ texture and store in _floatTextures_
        if (graphicsState.floatTextures->find(name) !=
            graphicsState.floatTextures->end())
            Warning("Texture \"%s\" being redefined", name.c_str());
        WARN_IF_ANIMATED_TRANSFORM("Texture");
        std::shared_ptr<Texture<Float>> ft =
            MakeFloatTexture(texname, curTransform[0], tp);
        if (ft) {
            // TODO: move this to be a GraphicsState method, also don't
            // provide direct floatTextures access?
            if (graphicsState.floatTexturesShared) {
                graphicsState.floatTextures =
                    std::make_shared<GraphicsState::FloatTextureMap>(*graphicsState.floatTextures);
                graphicsState.floatTexturesShared = false;
            }
            (*graphicsState.floatTextures)[name] = ft;
        }
    } else if (type == "color" || type == "spectrum") {
        // Create _color_ texture and store in _spectrumTextures_
        if (graphicsState.spectrumTextures->find(name) !=
            graphicsState.spectrumTextures->end())
            Warning("Texture \"%s\" being redefined", name.c_str());
        WARN_IF_ANIMATED_TRANSFORM("Texture");
        std::shared_ptr<Texture<Spectrum>> st =
            MakeSpectrumTexture(texname, curTransform[0], tp);
        if (st) {
            if (graphicsState.spectrumTexturesShared) {
                graphicsState.spectrumTextures =
                    std::make_shared<GraphicsState::SpectrumTextureMap>(*graphicsState.spectrumTextures);
                graphicsState.spectrumTexturesShared = false;
            }
            (*graphicsState.spectrumTextures)[name] = st;
        }
    } else
        Error("Texture type \"%s\" unknown.", type.c_str());
}

// 设置当前材质。通过 MakeMaterial 创建材质实例，并更新当前图形状态的材质信息。
void pbrtMaterial(const std::string &name, const ParamSet &params) {
    VERIFY_WORLD("Material");
    ParamSet emptyParams;
    TextureParams mp(params, emptyParams, *graphicsState.floatTextures,
                     *graphicsState.spectrumTextures);
    std::shared_ptr<Material> mtl = MakeMaterial(name, mp);
    graphicsState.currentMaterial =
        std::make_shared<MaterialInstance>(name, mtl, params);

    if (PbrtOptions.cat || PbrtOptions.toPly) {
        printf("%*sMaterial \"%s\" ", catIndentCount, "", name.c_str());
        params.Print(catIndentCount);
        printf("\n");
    }
}

// 创建命名材质，供后续 NamedMaterial 引用。材质类型在参数集的 "type" 字段中指定。
void pbrtMakeNamedMaterial(const std::string &name, const ParamSet &params) {
    VERIFY_WORLD("MakeNamedMaterial");
    // error checking, warning if replace, what to use for transform?
    ParamSet emptyParams;
    TextureParams mp(params, emptyParams, *graphicsState.floatTextures,
                     *graphicsState.spectrumTextures);
    std::string matName = mp.FindString("type");
    WARN_IF_ANIMATED_TRANSFORM("MakeNamedMaterial");
    if (matName == "")
        Error("No parameter string \"type\" found in MakeNamedMaterial");

    if (PbrtOptions.cat || PbrtOptions.toPly) {
        printf("%*sMakeNamedMaterial \"%s\" ", catIndentCount, "",
               name.c_str());
        params.Print(catIndentCount);
        printf("\n");
    } else {
        std::shared_ptr<Material> mtl = MakeMaterial(matName, mp);
        if (graphicsState.namedMaterials->find(name) !=
            graphicsState.namedMaterials->end())
            Warning("Named material \"%s\" redefined.", name.c_str());
        if (graphicsState.namedMaterialsShared) {
            graphicsState.namedMaterials =
                std::make_shared<GraphicsState::NamedMaterialMap>(*graphicsState.namedMaterials);
            graphicsState.namedMaterialsShared = false;
        }
        (*graphicsState.namedMaterials)[name] =
            std::make_shared<MaterialInstance>(matName, mtl, params);
    }
}

// 引用之前通过 MakeNamedMaterial 创建的命名材质，设置为当前材质。
void pbrtNamedMaterial(const std::string &name) {
    VERIFY_WORLD("NamedMaterial");
    if (PbrtOptions.cat || PbrtOptions.toPly) {
        printf("%*sNamedMaterial \"%s\"\n", catIndentCount, "", name.c_str());
        return;
    }

    auto iter = graphicsState.namedMaterials->find(name);
    if (iter == graphicsState.namedMaterials->end()) {
        Error("NamedMaterial \"%s\" unknown.", name.c_str());
        return;
    }
    graphicsState.currentMaterial = iter->second;
}

// 创建光源（非面光源）。通过 MakeLight 创建光源对象并添加到场景灯光列表中。
void pbrtLightSource(const std::string &name, const ParamSet &params) {
    VERIFY_WORLD("LightSource");
    WARN_IF_ANIMATED_TRANSFORM("LightSource");
    MediumInterface mi = graphicsState.CreateMediumInterface();
    std::shared_ptr<Light> lt = MakeLight(name, params, curTransform[0], mi);
    if (!lt)
        Error("LightSource: light type \"%s\" unknown.", name.c_str());
    else
        renderOptions->lights.push_back(lt);
    if (PbrtOptions.cat || PbrtOptions.toPly) {
        printf("%*sLightSource \"%s\" ", catIndentCount, "", name.c_str());
        params.Print(catIndentCount);
        printf("\n");
    }
}

// 设置面光源。面光源会与后续的 Shape 命令关联，使形状具有自发光属性。
void pbrtAreaLightSource(const std::string &name, const ParamSet &params) {
    VERIFY_WORLD("AreaLightSource");
    graphicsState.areaLight = name;
    graphicsState.areaLightParams = params;
    if (PbrtOptions.cat || PbrtOptions.toPly) {
        printf("%*sAreaLightSource \"%s\" ", catIndentCount, "", name.c_str());
        params.Print(catIndentCount);
        printf("\n");
    }
}

// 创建几何形状。处理静态和动画两种变换模式下的形状创建，
// 为形状分配材质，可选地创建面光源，并添加到场景图或当前实例中。
void pbrtShape(const std::string &name, const ParamSet &params) {
    VERIFY_WORLD("Shape");
    std::vector<std::shared_ptr<Primitive>> prims;
    std::vector<std::shared_ptr<AreaLight>> areaLights;
    if (PbrtOptions.cat || (PbrtOptions.toPly && name != "trianglemesh")) {
        printf("%*sShape \"%s\" ", catIndentCount, "", name.c_str());
        params.Print(catIndentCount);
        printf("\n");
    }

    if (!curTransform.IsAnimated()) {
        // Initialize _prims_ and _areaLights_ for static shape

        // Create shapes for shape _name_
        Transform *ObjToWorld = transformCache.Lookup(curTransform[0]);
        Transform *WorldToObj = transformCache.Lookup(Inverse(curTransform[0]));
        std::vector<std::shared_ptr<Shape>> shapes =
            MakeShapes(name, ObjToWorld, WorldToObj,
                       graphicsState.reverseOrientation, params);
        if (shapes.empty()) return;
        std::shared_ptr<Material> mtl = graphicsState.GetMaterialForShape(params);
        params.ReportUnused();
        MediumInterface mi = graphicsState.CreateMediumInterface();
        prims.reserve(shapes.size());
        for (auto s : shapes) {
            // Possibly create area light for shape
            std::shared_ptr<AreaLight> area;
            if (graphicsState.areaLight != "") {
                area = MakeAreaLight(graphicsState.areaLight, curTransform[0],
                                     mi, graphicsState.areaLightParams, s);
                if (area) areaLights.push_back(area);
            }
            prims.push_back(
                std::make_shared<GeometricPrimitive>(s, mtl, area, mi));
        }
    } else {
        // Initialize _prims_ and _areaLights_ for animated shape

        // Create initial shape or shapes for animated shape
        if (graphicsState.areaLight != "")
            Warning(
                "Ignoring currently set area light when creating "
                "animated shape");
        Transform *identity = transformCache.Lookup(Transform());
        std::vector<std::shared_ptr<Shape>> shapes = MakeShapes(
            name, identity, identity, graphicsState.reverseOrientation, params);
        if (shapes.empty()) return;

        // Create _GeometricPrimitive_(s) for animated shape
        std::shared_ptr<Material> mtl = graphicsState.GetMaterialForShape(params);
        params.ReportUnused();
        MediumInterface mi = graphicsState.CreateMediumInterface();
        prims.reserve(shapes.size());
        for (auto s : shapes)
            prims.push_back(
                std::make_shared<GeometricPrimitive>(s, mtl, nullptr, mi));

        // Create single _TransformedPrimitive_ for _prims_

        // Get _animatedObjectToWorld_ transform for shape
        static_assert(MaxTransforms == 2,
                      "TransformCache assumes only two transforms");
        Transform *ObjToWorld[2] = {
            transformCache.Lookup(curTransform[0]),
            transformCache.Lookup(curTransform[1])
        };
        AnimatedTransform animatedObjectToWorld(
            ObjToWorld[0], renderOptions->transformStartTime, ObjToWorld[1],
            renderOptions->transformEndTime);
        if (prims.size() > 1) {
            std::shared_ptr<Primitive> bvh = std::make_shared<BVHAccel>(prims);
            prims.clear();
            prims.push_back(bvh);
        }
        prims[0] = std::make_shared<TransformedPrimitive>(
            prims[0], animatedObjectToWorld);
    }
    // Add _prims_ and _areaLights_ to scene or current instance
    if (renderOptions->currentInstance) {
        if (areaLights.size())
            Warning("Area lights not supported with object instancing");
        renderOptions->currentInstance->insert(
            renderOptions->currentInstance->end(), prims.begin(), prims.end());
    } else {
        renderOptions->primitives.insert(renderOptions->primitives.end(),
                                         prims.begin(), prims.end());
        if (areaLights.size())
            renderOptions->lights.insert(renderOptions->lights.end(),
                                         areaLights.begin(), areaLights.end());
    }
}

// 启发式函数：判断形状的 ParamSet 是否可能包含材质参数。
// 由于材质没有显式参数描述，我们通过检查参数名称和数量来推断。
// 这用于决定是否为形状创建唯一的材质实例。
bool shapeMaySetMaterialParameters(const ParamSet &ps) {
    for (const auto &param : ps.textures)
        // Any texture other than one for an alpha mask is almost certainly
        // for a Material (or is unused!).
        if (param->name != "alpha" && param->name != "shadowalpha")
            return true;

    // Special case spheres, which are the most common non-mesh primitive.
    for (const auto &param : ps.floats)
        if (param->nValues == 1 && param->name != "radius")
            return true;

    // Extra special case strings, since plymesh uses "filename", curve "type",
    // and loopsubdiv "scheme".
    for (const auto &param : ps.strings)
        if (param->nValues == 1 && param->name != "filename" &&
            param->name != "type" && param->name != "scheme")
            return true;

    // For all other parameter types, if there is a single value of the
    // parameter, assume it may be for the material. This should be valid
    // (if conservative), since no materials currently take array
    // parameters.
    for (const auto &param : ps.bools)
        if (param->nValues == 1)
            return true;
    for (const auto &param : ps.ints)
        if (param->nValues == 1)
            return true;
    for (const auto &param : ps.point2fs)
        if (param->nValues == 1)
            return true;
    for (const auto &param : ps.vector2fs)
        if (param->nValues == 1)
            return true;
    for (const auto &param : ps.point3fs)
        if (param->nValues == 1)
            return true;
    for (const auto &param : ps.vector3fs)
        if (param->nValues == 1)
            return true;
    for (const auto &param : ps.normals)
        if (param->nValues == 1)
            return true;
    for (const auto &param : ps.spectra)
        if (param->nValues == 1)
            return true;

    return false;
}

// 为形状获取适当的材质。如果形状参数提供了材质参数值，
// 则创建新的材质实例，否则复用当前材质。
std::shared_ptr<Material> GraphicsState::GetMaterialForShape(
    const ParamSet &shapeParams) {
    CHECK(currentMaterial);
    if (shapeMaySetMaterialParameters(shapeParams)) {
        // Only create a unique material for the shape if the shape's
        // parameters are (apparently) going to provide values for some of
        // the material parameters.
        TextureParams mp(shapeParams, currentMaterial->params, *floatTextures,
                         *spectrumTextures);
        return MakeMaterial(currentMaterial->name, mp);
    } else
        return currentMaterial->material;
}

// 根据当前 GraphicsState 中设置的介质名称创建 MediumInterface。
// 在 renderOptions 的 namedMedia 表中查找对应的介质对象。
MediumInterface GraphicsState::CreateMediumInterface() {
    MediumInterface m;
    if (currentInsideMedium != "") {
        if (renderOptions->namedMedia.find(currentInsideMedium) !=
            renderOptions->namedMedia.end())
            m.inside = renderOptions->namedMedia[currentInsideMedium].get();
        else
            Error("Named medium \"%s\" undefined.",
                  currentInsideMedium.c_str());
    }
    if (currentOutsideMedium != "") {
        if (renderOptions->namedMedia.find(currentOutsideMedium) !=
            renderOptions->namedMedia.end())
            m.outside = renderOptions->namedMedia[currentOutsideMedium].get();
        else
            Error("Named medium \"%s\" undefined.",
                  currentOutsideMedium.c_str());
    }
    return m;
}

// 翻转法线朝向。切换当前图形状态的 reverseOrientation 标志，
// 用于修正几何体法线方向。
void pbrtReverseOrientation() {
    VERIFY_WORLD("ReverseOrientation");
    graphicsState.reverseOrientation = !graphicsState.reverseOrientation;
    if (PbrtOptions.cat || PbrtOptions.toPly)
        printf("%*sReverseOrientation\n", catIndentCount, "");
}

// 开始对象实例定义。保存当前属性状态，并在 renderOptions 中创建
// 当前实例的 Primitive 列表。
void pbrtObjectBegin(const std::string &name) {
    VERIFY_WORLD("ObjectBegin");
    pbrtAttributeBegin();
    if (renderOptions->currentInstance)
        Error("ObjectBegin called inside of instance definition");
    renderOptions->instances[name] = std::vector<std::shared_ptr<Primitive>>();
    renderOptions->currentInstance = &renderOptions->instances[name];
    if (PbrtOptions.cat || PbrtOptions.toPly)
        printf("%*sObjectBegin \"%s\"\n", catIndentCount, "", name.c_str());
}

STAT_COUNTER("Scene/Object instances created", nObjectInstancesCreated);

// 结束对象实例定义。恢复之前保存的属性状态。
void pbrtObjectEnd() {
    VERIFY_WORLD("ObjectEnd");
    if (!renderOptions->currentInstance)
        Error("ObjectEnd called outside of instance definition");
    if (PbrtOptions.cat || PbrtOptions.toPly)
        printf("%*sObjectEnd\n", catIndentCount, "");
    renderOptions->currentInstance = nullptr;
    pbrtAttributeEnd();
    ++nObjectInstancesCreated;
}

STAT_COUNTER("Scene/Object instances used", nObjectInstancesUsed);

// 在当前变换下实例化之前通过 ObjectBegin/ObjectEnd 定义的对象。
// 为实例创建 TransformedPrimitive 包装器。
void pbrtObjectInstance(const std::string &name) {
    VERIFY_WORLD("ObjectInstance");
    if (PbrtOptions.cat || PbrtOptions.toPly) {
        printf("%*sObjectInstance \"%s\"\n", catIndentCount, "", name.c_str());
        return;
    }

    // Perform object instance error checking
    if (renderOptions->currentInstance) {
        Error("ObjectInstance can't be called inside instance definition");
        return;
    }
    if (renderOptions->instances.find(name) == renderOptions->instances.end()) {
        Error("Unable to find instance named \"%s\"", name.c_str());
        return;
    }
    std::vector<std::shared_ptr<Primitive>> &in =
        renderOptions->instances[name];
    if (in.empty()) return;
    ++nObjectInstancesUsed;
    if (in.size() > 1) {
        // Create aggregate for instance _Primitive_s
        std::shared_ptr<Primitive> accel(
            MakeAccelerator(renderOptions->AcceleratorName, std::move(in),
                            renderOptions->AcceleratorParams));
        if (!accel) accel = std::make_shared<BVHAccel>(in);
        in.clear();
        in.push_back(accel);
    }
    static_assert(MaxTransforms == 2,
                  "TransformCache assumes only two transforms");
    // Create _animatedInstanceToWorld_ transform for instance
    Transform *InstanceToWorld[2] = {
        transformCache.Lookup(curTransform[0]),
        transformCache.Lookup(curTransform[1])
    };
    AnimatedTransform animatedInstanceToWorld(
        InstanceToWorld[0], renderOptions->transformStartTime,
        InstanceToWorld[1], renderOptions->transformEndTime);
    std::shared_ptr<Primitive> prim(
        std::make_shared<TransformedPrimitive>(in[0], animatedInstanceToWorld));
    renderOptions->primitives.push_back(prim);
}

// 结束世界块。创建积分器和场景，调用积分器的 Render 方法进行渲染。
// 渲染完成后清理资源、报告统计信息。
void pbrtWorldEnd() {
    VERIFY_WORLD("WorldEnd");
    // Ensure there are no pushed graphics states
    while (pushedGraphicsStates.size()) {
        Warning("Missing end to pbrtAttributeBegin()");
        pushedGraphicsStates.pop_back();
        pushedTransforms.pop_back();
    }
    while (pushedTransforms.size()) {
        Warning("Missing end to pbrtTransformBegin()");
        pushedTransforms.pop_back();
    }

    // Create scene and render
    if (PbrtOptions.cat || PbrtOptions.toPly) {
        printf("%*sWorldEnd\n", catIndentCount, "");
    } else {
        std::unique_ptr<Integrator> integrator(renderOptions->MakeIntegrator());
        std::unique_ptr<Scene> scene(renderOptions->MakeScene());

        // This is kind of ugly; we directly override the current profiler
        // state to switch from parsing/scene construction related stuff to
        // rendering stuff and then switch it back below. The underlying
        // issue is that all the rest of the profiling system assumes
        // hierarchical inheritance of profiling state; this is the only
        // place where that isn't the case.
        CHECK_EQ(CurrentProfilerState(), ProfToBits(Prof::SceneConstruction));
        ProfilerState = ProfToBits(Prof::IntegratorRender);

        if (scene && integrator) integrator->Render(*scene);

        CHECK_EQ(CurrentProfilerState(), ProfToBits(Prof::IntegratorRender));
        ProfilerState = ProfToBits(Prof::SceneConstruction);
    }

    // Clean up after rendering. Do this before reporting stats so that
    // destructors can run and update stats as needed.
    graphicsState = GraphicsState();
    transformCache.Clear();
    currentApiState = APIState::OptionsBlock;
    ImageTexture<Float, Float>::ClearCache();
    ImageTexture<RGBSpectrum, Spectrum>::ClearCache();
    renderOptions.reset(new RenderOptions);

    if (!PbrtOptions.cat && !PbrtOptions.toPly) {
        MergeWorkerThreadStats();
        ReportThreadStats();
        if (!PbrtOptions.quiet) {
            PrintStats(stdout);
            ReportProfilerResults(stdout);
            ClearStats();
            ClearProfiler();
        }
    }

    for (int i = 0; i < MaxTransforms; ++i) curTransform[i] = Transform();
    activeTransformBits = AllTransformsBits;
    namedCoordinateSystems.erase(namedCoordinateSystems.begin(),
                                 namedCoordinateSystems.end());
}

// 创建场景对象。使用配置好的加速器将所有图元组织为加速结构，
// 并结合所有已定义的光源构造 Scene。
Scene *RenderOptions::MakeScene() {
    std::shared_ptr<Primitive> accelerator =
        MakeAccelerator(AcceleratorName, std::move(primitives), AcceleratorParams);
    if (!accelerator) accelerator = std::make_shared<BVHAccel>(primitives);
    Scene *scene = new Scene(accelerator, lights);
    // Erase primitives and lights from _RenderOptions_
    primitives.clear();
    lights.clear();
    return scene;
}

// 创建积分器对象。先创建相机和采样器，然后根据配置的积分器名称
// 创建对应的积分器实例（whitted、directlighting、path、volpath、bdpt、mlt、ao、sppm 等）。
Integrator *RenderOptions::MakeIntegrator() const {
    std::shared_ptr<const Camera> camera(MakeCamera());
    if (!camera) {
        Error("Unable to create camera");
        return nullptr;
    }

    std::shared_ptr<Sampler> sampler =
        MakeSampler(SamplerName, SamplerParams, camera->film);
    if (!sampler) {
        Error("Unable to create sampler.");
        return nullptr;
    }

    Integrator *integrator = nullptr;
    if (IntegratorName == "whitted")
        integrator = CreateWhittedIntegrator(IntegratorParams, sampler, camera);
    else if (IntegratorName == "directlighting")
        integrator =
            CreateDirectLightingIntegrator(IntegratorParams, sampler, camera);
    else if (IntegratorName == "path")
        integrator = CreatePathIntegrator(IntegratorParams, sampler, camera);
    else if (IntegratorName == "volpath")
        integrator = CreateVolPathIntegrator(IntegratorParams, sampler, camera);
    else if (IntegratorName == "bdpt") {
        integrator = CreateBDPTIntegrator(IntegratorParams, sampler, camera);
    } else if (IntegratorName == "mlt") {
        integrator = CreateMLTIntegrator(IntegratorParams, camera);
    } else if (IntegratorName == "ambientocclusion") {
        integrator = CreateAOIntegrator(IntegratorParams, sampler, camera);
    } else if (IntegratorName == "sppm") {
        integrator = CreateSPPMIntegrator(IntegratorParams, camera);
    } else {
        Error("Integrator \"%s\" unknown.", IntegratorName.c_str());
        return nullptr;
    }

    if (renderOptions->haveScatteringMedia && IntegratorName != "volpath" &&
        IntegratorName != "bdpt" && IntegratorName != "mlt") {
        Warning(
            "Scene has scattering media but \"%s\" integrator doesn't support "
            "volume scattering. Consider using \"volpath\", \"bdpt\", or "
            "\"mlt\".", IntegratorName.c_str());
    }

    IntegratorParams.ReportUnused();
    // Warn if no light sources are defined
    if (lights.empty())
        Warning(
            "No light sources defined in scene; "
            "rendering a black image.");
    return integrator;
}

// 创建相机对象。依次创建滤波器、胶片，然后调用 MakeCamera 构造相机。
Camera *RenderOptions::MakeCamera() const {
    std::unique_ptr<Filter> filter = MakeFilter(FilterName, FilterParams);
    Film *film = MakeFilm(FilmName, FilmParams, std::move(filter));
    if (!film) {
        Error("Unable to create film.");
        return nullptr;
    }
    Camera *camera = pbrt::MakeCamera(CameraName, CameraParams, CameraToWorld,
                                  renderOptions->transformStartTime,
                                  renderOptions->transformEndTime, film);
    return camera;
}

}  // namespace pbrt
