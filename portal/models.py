from django.db import models


class Question(models.Model):
    SUBJECT_CHOICES = [
        ("physics", "Physics"),
        ("maths", "Maths"),
        ("chemistry", "Chemistry"),
    ]

    subject = models.CharField(max_length=20, choices=SUBJECT_CHOICES)
    topic = models.CharField(max_length=100, blank=True)
    text = models.TextField()
    created_at = models.DateTimeField(auto_now_add=True)
    status = models.CharField(
        max_length=20,
        default="new",  # later: new / in-progress / solved
    )

    def __str__(self) -> str:
        return f"{self.subject} – {self.topic or 'No topic'}"
