import spiceypy as spice

spice.furnsh("naif0012.tls")
print("TLS loaded")

spice.furnsh("de440s.bsp")
print("BSP loaded")

print("SPICE working")