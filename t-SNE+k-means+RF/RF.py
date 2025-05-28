import pandas as pd
from sklearn.model_selection import train_test_split, StratifiedKFold, GridSearchCV, cross_val_score
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import confusion_matrix, classification_report
from sklearn.manifold import TSNE
from sklearn.cluster import KMeans
import matplotlib
matplotlib.use('Agg')  # 选择非交互式后端
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import os
from joblib import parallel_backend
# 读取并转置数据
df = pd.read_csv("D:\\Bioana\\dbparsemergedv1\\countDatamerged2.csv",
                 sep='\t', index_col=0, encoding='utf-8').transpose()

# 分配类别标签：crispatus -> 1, iners -> 0
df['Lactotype'] = df.index.str.contains('crispatus').astype(int)
# 创建一个临时数据框，将 'Lactotype' 移动到第一列
cols = ['Lactotype'] + [col for col in df.columns if col != 'Lactotype']
dftemp = df[cols]

# 创建输出目录并保存临时数据框到 CSV 文件
save_path = 'new/df.csv'
os.makedirs(os.path.dirname(save_path), exist_ok=True)  # 如果目录不存在则创建
dftemp.to_csv(save_path, sep='\t', encoding='utf-8')

print(f"分类后的数据已保存到 {save_path}")
# 特征和标签
X = df.drop('Lactotype', axis=1).values
y = df['Lactotype'].values

# 拆分数据集
X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.3, random_state=42, stratify=y)

# 使用 t-SNE 降维
tsne = TSNE(n_components=2, learning_rate='auto', init='random', perplexity=30, random_state=42)
X_tsne = tsne.fit_transform(X_train)

# 使用 KMeans 聚类
kmeans = KMeans(n_clusters=2, random_state=42)
cluster_labels = kmeans.fit_predict(X_train)

# 可视化 t-SNE 和 KMeans 聚类结果
plt.figure(figsize=(12, 8))
scatter = plt.scatter(X_tsne[:, 0], X_tsne[:, 1], c=cluster_labels, cmap='viridis', s=30)
plt.colorbar(scatter, label='Cluster Label')
plt.title('t-SNE Visualization with KMeans Clustering')
plt.xlabel('t-SNE Component 1')
plt.ylabel('t-SNE Component 2')
save_path = 'new/tsne_kmeans.png'
os.makedirs(os.path.dirname(save_path), exist_ok=True)  # 如果目录不存在则创建
plt.savefig(save_path)
plt.close()

# 定义随机森林分类器
rfc = RandomForestClassifier(random_state=42, class_weight='balanced')

# 设置GridSearchCV的超参数范围
param_grid = {
    'n_estimators': [100],
    'max_depth': [8],
    'min_samples_split': [20],
    'min_samples_leaf': [20],
    'max_features': ['sqrt']
}
# 在代码开头创建一个只包含 ASCII 路径的临时文件夹
TEMP_FOLDER = r"C:\temp_joblib"
os.makedirs(TEMP_FOLDER, exist_ok=True)

# 在使用 GridSearchCV 前通过 parallel_backend 显式指定 temp_folder
with parallel_backend('loky', n_jobs=-1, temp_folder=TEMP_FOLDER):
    grid_search = GridSearchCV(
        estimator=rfc,
        param_grid=param_grid,
        cv=5,
        scoring='accuracy',
        # 注意：这里不再写 n_jobs=-1 ，因为我们在 parallel_backend 已经指定了并行
        verbose=2
    )
    grid_search.fit(X_train, y_train)

# 输出最佳参数
print(f"Best parameters from GridSearchCV: {grid_search.best_params_}")

# 使用最佳参数训练模型
best_rfc = grid_search.best_estimator_

# 十折交叉验证
cv = StratifiedKFold(n_splits=10, shuffle=True, random_state=42)
scores = cross_val_score(best_rfc, X_train, y_train, cv=cv, scoring='accuracy')

print(f"Accuracy for each fold: {scores}")
print(f"Mean accuracy across all folds: {scores.mean():.4f}")

# 使用测试集评估模型
best_rfc.fit(X_train, y_train)
prediction = best_rfc.predict(X_test)
print(f"Test Accuracy: {best_rfc.score(X_test, y_test):.4f}")

# 输出分类报告
print(f"Classification report:\n{classification_report(y_test, prediction)}")

# 绘制混淆矩阵
def plot_confusion_matrix(y_true, y_pred):
    fig, ax = plt.subplots(figsize=(8, 6))
    cm = confusion_matrix(y_true, y_pred, labels=[1, 0])
    print(cm)
    # 计算混淆矩阵中的四个元素
    TP = cm[0, 0]  # True Positive (True Crispatus)
    FN = cm[0, 1]  # False Negative (False Crispatus)
    FP = cm[1, 0]  # False Positive (False Iners)
    TN = cm[1, 1]  # True Negative (True Iners)
    print(f"FP: {FP}")
    # 计算各项指标
    TPR = TP / (TP + FN) if (TP + FN) > 0 else 0  # 真阳性率
    FNR = FN / (TP + FN) if (TP + FN) > 0 else 0  # 假阴性率
    FPR = FP / (FP + TN) if (FP + TN) > 0 else 0  # 假阳性率
    TNR = TN / (FP + TN) if (FP + TN) > 0 else 0  # 真阴性率
    # 更新混淆矩阵标签
    group_names = ["True Crispatus", "False Iners", "False Crispatus", "True Iners"]
    group_counts = ["{0:0.0f}".format(value) for value in cm.flatten()]
    # 更新group_percentages，使用这些指标来代替原先的百分比
    group_percentages = [
        "{0:.2%}".format(TPR),  # True Crispatus
        "{0:.2%}".format(FNR),  # False Iners
        "{0:.2%}".format(FPR),  # False Crispatus
        "{0:.2%}".format(TNR)  # True Iners
    ]

    labels = [f"{v1}\n{v2}\n{v3}" for v1, v2, v3 in zip(group_names, group_counts, group_percentages)]
    labels = np.asarray(labels).reshape(2, 2)

    # 绘制热力图
    sns.set(font_scale=1.4, style='whitegrid', palette='pastel')
    sns.heatmap(cm, annot=labels, fmt="", cmap="Blues_r", cbar=False, ax=ax)

    ax.set_xlabel('Predicted Class', fontsize=16)
    ax.set_ylabel('Actual Class', fontsize=16)
    ax.set_xticklabels(['Crispatus', 'Iners'], fontsize=14)
    ax.set_yticklabels(['Crispatus', 'Iners'], fontsize=14)
    ax.set_title('Confusion Matrix - Random Forest Classification', fontsize=16)

    # 创建目录 new 并保存图片
    save_path = 'new/confusion_matrix.png'
    os.makedirs(os.path.dirname(save_path), exist_ok=True)  # 如果目录不存在则创建
    plt.savefig(save_path)
    plt.close()


plot_confusion_matrix(y_test, prediction)

# 特征重要性
importances = best_rfc.feature_importances_
indices = np.argsort(importances)[::-1]

# 输出前一半的特征重要性
num_features_to_display = len(importances) // 5  # 前1/5
top_features = indices[:num_features_to_display]

# 可视化前一半特征的重要性
plt.figure(figsize=(14, 8))
plt.title("Top Feature Importances")
plt.bar(range(num_features_to_display), importances[top_features], align="center")
plt.xticks(range(num_features_to_display), df.columns[top_features], rotation=90)
plt.xlabel('Features')
plt.ylabel('Importance Score')
plt.tight_layout()
save_path = 'new/top_feature_importance.png'
os.makedirs(os.path.dirname(save_path), exist_ok=True)  # 如果目录不存在则创建
plt.savefig(save_path)
plt.close()
