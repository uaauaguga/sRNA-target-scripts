from torch import nn
from torch.nn import functional as F
import torch
import numpy as np

class ResidualBlock2D(nn.Module):
    def __init__(self, planes, kernel_size=(5,5), padding=(2,2), downsample=True):
        super(ResidualBlock2D, self).__init__()
        self.conv1 = nn.Conv2d(planes,   planes,   kernel_size=1, stride=1, bias=False)
        self.bn1 = nn.BatchNorm2d(planes)
        self.conv2 = nn.Conv2d(planes,   planes*2, kernel_size=kernel_size, stride=1,
                     padding=padding, bias=False)
        self.bn2 = nn.BatchNorm2d(planes*2)
        self.conv3 = nn.Conv2d(planes*2, planes*4, kernel_size=1, stride=1, bias=False)
        self.bn3 = nn.BatchNorm2d(planes * 4)
        self.downsample = nn.Sequential(
            nn.Conv2d(planes,   planes*4,   kernel_size=1, stride=1, bias=False),
            nn.BatchNorm2d(planes*4),
        )
        self.relu  = nn.ReLU(inplace=True)

    def forward(self, x):
        identity = x

        out = self.conv1(x)
        out = self.bn1(out)
        out = self.relu(out)

        out = self.conv2(out)
        out = self.bn2(out)
        out = self.relu(out)

        out = self.conv3(out)
        out = self.bn3(out)

        if self.downsample:
            identity = self.downsample(x)
        out += identity
        out = self.relu(out)

        return out

class InteractionScorer(nn.Module):
    def __init__(self,planes=8):
        super(InteractionScorer, self).__init__()
        self.planes = planes
        self.conv = nn.Conv2d(8, self.planes, kernel_size=7, stride=2, padding=3, bias=False)
        self.bn = nn.BatchNorm2d(planes)
        self.maxpool = nn.MaxPool2d(kernel_size=3, stride=2, padding=1)
        self.res1 = ResidualBlock2D(self.planes)
        self.res2 = ResidualBlock2D(self.planes*4)
        self.res3 = ResidualBlock2D(self.planes*16)
        self.avgpool = nn.AdaptiveAvgPool2d((1, 1))
        self.fc = nn.Linear(self.planes*64, 1)
    def forward(self,x):
        x = self.conv(x)
        x = self.bn(x).relu()
        x = self.maxpool(x)
        x = self.res1(x)
        x = self.res2(x)
        x = self.res3(x)
        x = self.avgpool(x).flatten(1)
        x = self.fc(x).sigmoid()
        return x

