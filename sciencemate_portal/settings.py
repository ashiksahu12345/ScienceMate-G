"""
Django settings for sciencemate_portal project.
"""

from pathlib import Path
import os

# ---------------- Base Paths ----------------
BASE_DIR = Path(__file__).resolve().parent.parent

# ---------------- Security ------------------
SECRET_KEY = "django-insecure-123456replace-this-when-production"

DEBUG = True   # Render पर बाद में False भी कर सकते हो
ALLOWED_HOSTS = ["*"]   # अभी के लिए सब allow

# ---------------- Installed Apps ------------
INSTALLED_APPS = [
    "django.contrib.admin",
    "django.contrib.auth",
    "django.contrib.contenttypes",
    "django.contrib.sessions",
    "django.contrib.messages",
    "django.contrib.staticfiles",

    # main app
    "portal",
]

# ---------------- Middleware ----------------
MIDDLEWARE = [
    "django.middleware.security.SecurityMiddleware",
    "whitenoise.middleware.WhiteNoiseMiddleware",       # static files for Render
    "django.contrib.sessions.middleware.SessionMiddleware",
    "django.middleware.common.CommonMiddleware",
    "django.middleware.csrf.CsrfViewMiddleware",
    "django.contrib.auth.middleware.AuthenticationMiddleware",
    "django.contrib.messages.middleware.MessageMiddleware",
    "django.middleware.clickjacking.XFrameOptionsMiddleware",
]

# ---------------- URLs / WSGI ---------------
ROOT_URLCONF = "sciencemate_portal.urls"
WSGI_APPLICATION = "sciencemate_portal.wsgi.application"

# ---------------- Templates -----------------
TEMPLATES = [
    {
        "BACKEND": "django.template.backends.django.DjangoTemplates",
        "DIRS": [
            BASE_DIR / "templates",   # global templates folder
        ],
        "APP_DIRS": True,
        "OPTIONS": {
            "context_processors": [
                "django.template.context_processors.debug",
                "django.template.context_processors.request",
                "django.contrib.auth.context_processors.auth",
                "django.contrib.messages.context_processors.messages",
            ],
        },
    },
]

# ---------------- Database ------------------
DATABASES = {
    "default": {
        "ENGINE": "django.db.backends.sqlite3",
        "NAME": BASE_DIR / "db.sqlite3",
    }
}

# -------- Password Validation ---------------
AUTH_PASSWORD_VALIDATORS = [
    {"NAME": "django.contrib.auth.password_validation.UserAttributeSimilarityValidator"},
    {"NAME": "django.contrib.auth.password_validation.MinimumLengthValidator"},
    {"NAME": "django.contrib.auth.password_validation.CommonPasswordValidator"},
    {"NAME": "django.contrib.auth.password_validation.NumericPasswordValidator"},
]

# --------- i18n / Timezone ------------------
LANGUAGE_CODE = "en-us"
TIME_ZONE = "Asia/Kolkata"
USE_I18N = True
USE_TZ = True

# ---------------- Static Files --------------
STATIC_URL = "/static/"

# यहाँ वो folder दो जहाँ तुम्हारे original CSS/JS रखे हैं
# आपके project में screenshot से "static_portal" दिख रहा था
STATICFILES_DIRS = [
    BASE_DIR / "static",
]

# collectstatic के बाद Render यहीं से serve करेगा
STATIC_ROOT = BASE_DIR / "staticfiles"

# Whitenoise storage (compressed + hashed filenames)
STATICFILES_STORAGE = "whitenoise.storage.CompressedManifestStaticFilesStorage"

# ---------------- Media Files ---------------
MEDIA_URL = "/media/"
MEDIA_ROOT = BASE_DIR / "media"

# ------------- Default Primary Key ----------
DEFAULT_AUTO_FIELD = "django.db.models.BigAutoField"
