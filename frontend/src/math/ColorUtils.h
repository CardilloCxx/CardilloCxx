namespace frontend::math {
    inline filament::math::float4 encodeEntityPickColor(uint32_t entityId) {
        const uint32_t encoded = entityId + 1u;
        return {
            static_cast<float>((encoded >> 16) & 0xff) / 255.0f,
            static_cast<float>((encoded >> 8) & 0xff) / 255.0f,
            static_cast<float>(encoded & 0xff) / 255.0f,
            1.0f
        };
    }
}