import pandas as pd
import plotly.express as px
import numpy as np
import plotly.colors
import os

# 读取CSV数据
df = pd.read_csv('1CSE.csv')  # 修改为实际文件路径

# 解析突变信息
positions = []  # 存储所有位置标签（保持顺序）
data = []
for _, row in df.iterrows():
    mutation = row['mutation']
    ddg = row['ddG_pred']
    
    # 解析字段：例如 "A_E_0_C" -> wt=A, chain=E, pos=0, mut=C
    parts = mutation.split('_')
    wt = parts[0]
    chain = parts[1]
    pos = parts[2]
    mut = parts[3]
    
    # 跳过 C、G、P 的突变
    if mut in ['C', 'G', 'P']:
        continue
    
    # 构建位置标签（例如 "A0"）
    pos_label = f"{wt}{pos}_{chain}"
    
    # 记录位置顺序（确保唯一性）
    if pos_label not in positions:
        positions.append(pos_label)
    
    data.append([pos_label, mut, ddg])

# 转换为DataFrame
df_processed = pd.DataFrame(data, columns=['Position', 'Mutation', 'ddG'])

# 定义要显示的位置区间（根据位置数字过滤）
selected_positions = []
for pos_label in positions:
    # 提取位置数字（例如从 "A0" 提取 0）
    pos_num = int(''.join(filter(str.isdigit, pos_label)))
    
    # 筛选区间：28-35, 49-67, 100-111（根据实际需求调整）
    if 1<0:#(28 <= pos_num <= 35) or (49 <= pos_num <= 67) or (100 <= pos_num <= 111):
        selected_positions.append(pos_label)
    else:
        selected_positions.append(pos_label)

# 过滤数据，只保留选定的位置
df_filtered = df_processed[df_processed['Position'].isin(selected_positions)]

# 创建矩阵（行：突变氨基酸，列：位置，值：ddG）
amino_acids = sorted(df_filtered['Mutation'].unique())  # 突变氨基酸按字母排序
matrix = pd.DataFrame(index=amino_acids, columns=selected_positions)

# 填充矩阵，保留2位小数
for _, row in df_filtered.iterrows():
    matrix.at[row['Mutation'], row['Position']] = round(row['ddG'], 2)

# 获取数据范围（用于颜色比例）
zmin = np.nanmin(matrix.values.astype(float))
zmax = np.nanmax(matrix.values.astype(float))

# 自定义颜色比例（蓝-白-红）
zero_pos = (0 - zmin) / (zmax - zmin)  # 计算0点的位置
blue_gradient = plotly.colors.sample_colorscale('Blues', np.linspace(1, 0, 5))  # 深蓝到白
red_gradient = plotly.colors.sample_colorscale('Reds', np.linspace(0, 1, 5))    # 白到深红
custom_colorscale = []
for color in blue_gradient:
    custom_colorscale.append(color)
for color in red_gradient[1:]:  # 避免重复白色
    custom_colorscale.append(color)

# 生成热力图
fig = px.imshow(
    matrix.astype(float),
    labels=dict(x="Position", y="Mutated to", color="ΔΔG (kcal/mol)"),
    color_continuous_scale=custom_colorscale,
    aspect="auto",
    zmin=zmin,
    zmax=zmax,
    text_auto=True
)

# 更新布局：显示所有 x 轴标签并启用滚动条
fig.update_layout(
    title="Protein Mutation ΔΔG Predictions (1CSE)",
    width=2000,  # 增加宽度以适应长 x 轴
    height=600,
    xaxis=dict(
        tickmode='array',
        tickvals=list(range(len(selected_positions))),
        #ticktext=selected_positions,  # 显示所有标签
        ticktext=[selected_positions[i] if i % 3 == 0 else '' for i in range(len(selected_positions))],  # 每隔3个显示一个标签
        tickangle=45,  # 标签倾斜角度
        showgrid=False,
        fixedrange=False  # 允许用户缩放和拖动
    ),
    yaxis=dict(tickmode='array', tickvals=list(range(len(amino_acids)))),
    margin=dict(l=50, r=50, t=50, b=50),  # 调整边距
    template="plotly_white"
)

# 将图表保存为HTML，并添加滚动容器
output_path = os.path.join(os.getcwd(), '1CSE_mutation_ddG_heatmap.html')
html_content = f"""
<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>Heatmap with Scroll</title>
    <style>
        .scroll-container {{
            width: 1200px; /* 固定宽度 */
            overflow-x: auto; /* 启用水平滚动 */
            border: 1px solid #ccc;
            padding: 10px;
        }}
    </style>
</head>
<body>
    <h1>Protein Mutation ΔΔG Predictions (1CSE)</h1>
    <div class="scroll-container">
        {fig.to_html(full_html=False, include_plotlyjs='cdn')}
    </div>
</body>
</html>
"""

# 写入HTML文件
with open(output_path, 'w', encoding='utf-8') as f:
    f.write(html_content)

print(f"热力图已保存至：{output_path}")

# 显示图表
fig.show()