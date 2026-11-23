import numpy as np
import pandas as pd
from sklearn.model_selection import StratifiedKFold, LeaveOneOut
from sklearn.linear_model import LogisticRegression
from sklearn.ensemble import RandomForestClassifier
from sklearn.svm import SVC
from sklearn.preprocessing import StandardScaler
from sklearn.pipeline import Pipeline
from sklearn.metrics import (
    accuracy_score, precision_score, recall_score, f1_score,
    roc_auc_score, roc_curve, confusion_matrix
)
import plotly.graph_objects as go
import json

def get_classifier(classifier_type):
    """Factory to instantiate classifiers with standard hyperparameters."""
    clf_type = classifier_type.lower()
    if "forest" in clf_type or clf_type == "rf":
        return RandomForestClassifier(n_estimators=100, random_state=42, max_depth=5), "Random Forest", "#0d9488"
    elif "svm" in clf_type or "svc" in clf_type or "support_vector" in clf_type:
        return SVC(kernel='linear', probability=True, random_state=42, C=1.0), "Support Vector Machine (SVM)", "#06b6d4"
    else:
        # Default Logistic Regression (L2 penalty)
        return LogisticRegression(C=1.0, solver='liblinear', random_state=42), "Logistic Regression", "#2563eb"

def evaluate_classifier_cv(X, y, estimator, name, color, cv, feature_names):
    """
    Evaluates an estimator via cross-validation with strict data-leakage safeguards.
    The StandardScaler is strictly instantiated and fit inside each training fold!
    """
    y_true_all = []
    y_pred_all = []
    y_prob_all = []
    fold_accuracies = []

    for fold_idx, (train_idx, test_idx) in enumerate(cv.split(X, y)):
        X_train, X_test = X[train_idx], X[test_idx]
        y_train, y_test = y[train_idx], y[test_idx]

        # Scaler fitted ONLY on X_train to prevent leakage!
        pipeline = Pipeline([
            ('scaler', StandardScaler()),
            ('classifier', estimator)
        ])
        pipeline.fit(X_train, y_train)

        y_pred = pipeline.predict(X_test)
        if hasattr(pipeline.named_steps['classifier'], 'predict_proba'):
            y_prob = pipeline.predict_proba(X_test)[:, 1]
        elif hasattr(pipeline.named_steps['classifier'], 'decision_function'):
            y_prob = pipeline.decision_function(X_test)
        else:
            y_prob = y_pred.astype(float)

        y_true_all.extend(y_test)
        y_pred_all.extend(y_pred)
        y_prob_all.extend(y_prob)
        fold_accuracies.append(accuracy_score(y_test, y_pred))

    y_true_all = np.array(y_true_all)
    y_pred_all = np.array(y_pred_all)
    y_prob_all = np.array(y_prob_all)

    acc = float(accuracy_score(y_true_all, y_pred_all))
    prec = float(precision_score(y_true_all, y_pred_all, zero_division=0))
    rec = float(recall_score(y_true_all, y_pred_all, zero_division=0))
    f1 = float(f1_score(y_true_all, y_pred_all, zero_division=0))

    try:
        auc = float(roc_auc_score(y_true_all, y_prob_all))
        fpr_v, tpr_v, _ = roc_curve(y_true_all, y_prob_all)
    except Exception:
        auc = 0.5
        fpr_v, tpr_v = [0, 1], [0, 1]

    cm = confusion_matrix(y_true_all, y_pred_all)
    cm_dict = {
        'tn': int(cm[0, 0]) if cm.shape == (2, 2) else 0,
        'fp': int(cm[0, 1]) if cm.shape == (2, 2) else 0,
        'fn': int(cm[1, 0]) if cm.shape == (2, 2) else 0,
        'tp': int(cm[1, 1]) if cm.shape == (2, 2) else 0
    }

    # Fit final model on full training set for feature importances / independent cohort testing
    final_pipeline = Pipeline([
        ('scaler', StandardScaler()),
        ('classifier', estimator)
    ])
    final_pipeline.fit(X, y)

    gene_importances = []
    clf_obj = final_pipeline.named_steps['classifier']
    if hasattr(clf_obj, 'coef_'):
        coefs = clf_obj.coef_[0]
        for g, c in zip(feature_names, coefs):
            gene_importances.append({
                'Gene': g,
                'Importance': round(float(abs(c)), 4),
                'Direction': 'Positive' if c > 0 else 'Negative',
                'Raw_Weight': round(float(c), 4)
            })
    elif hasattr(clf_obj, 'feature_importances_'):
        imps = clf_obj.feature_importances_
        for g, i in zip(feature_names, imps):
            gene_importances.append({
                'Gene': g,
                'Importance': round(float(i), 4),
                'Direction': 'Nonlinear Importance',
                'Raw_Weight': round(float(i), 4)
            })

    gene_importances.sort(key=lambda x: x['Importance'], reverse=True)

    return {
        'model_name': name,
        'color': color,
        'accuracy': acc,
        'precision': prec,
        'recall': rec,
        'f1_score': f1,
        'roc_auc': auc,
        'fpr': list(fpr_v),
        'tpr': list(tpr_v),
        'confusion_matrix': cm_dict,
        'fold_accuracies': [round(a, 3) for a in fold_accuracies],
        'gene_importances': gene_importances,
        'pipeline': final_pipeline
    }
