from django import forms
from .models import Question


class ChatForm(forms.Form):
    message = forms.CharField(
        label="Your message",
        widget=forms.Textarea(
            attrs={
                "rows": 3,
                "class": "form-control",
                "placeholder": "Ask anything about Physics, Maths or Chemistry...",
            }
        ),
    )


class QuestionForm(forms.ModelForm):
    class Meta:
        model = Question
        fields = ["subject", "topic", "text"]
        widgets = {
            "subject": forms.Select(attrs={"class": "form-select"}),
            "topic": forms.TextInput(
                attrs={
                    "class": "form-control",
                    "placeholder": "e.g. Motion, Algebra, Atomic structure",
                }
            ),
            "text": forms.Textarea(
                attrs={
                    "class": "form-control",
                    "rows": 5,
                    "placeholder": "Paste or type your full question here...",
                }
            ),
        }
