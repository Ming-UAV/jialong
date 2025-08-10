

from PIL import Image

# 1. 读取底图和点图，并转换为 RGBA 模式
map_img = Image.open(r"C:\Users\28492\OneDrive\图片\Screenshots\Screenshot 2025-08-08 151139.png").convert("RGBA")
pts_img = Image.open(r"C:\Users\28492\OneDrive\Desktop\dissertation-GA-method\results\c104_20250706_123018.png").convert("RGBA")

# 2. 若两图尺寸不一致，可按底图尺寸缩放点图
pts_img = pts_img.resize(map_img.size, Image.LANCZOS)

# 3. 根据需要调整点图透明度（可选），此处设置为 80% 不透明
alpha = 100  # 0-255，255 完全不透明
r, g, b, a = pts_img.split()
pts_img.putalpha(Image.new("L", pts_img.size, alpha))

# 4. 叠加合成
combined = Image.alpha_composite(map_img, pts_img)

# 5. 保存与展示
import os

out_dir = r"D:\OneDrive\Desktop\dissertation-GA-method\results"
os.makedirs(out_dir, exist_ok=True)
out_path = os.path.join(out_dir, "overlay_serengeti_uav.png")

combined.convert("RGB").save(out_path, dpi=(300,300))

combined.show()