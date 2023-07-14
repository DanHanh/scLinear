"""Babel and neural network models implemented in PyTorch."""
import pytorch_lightning as pl
import torch
import torch.nn as nn
from torch.optim import Adam


class IdentityTransformer:
    @staticmethod
    def fit_transform(X):
        return X

    @staticmethod
    def transform(X):
        return X


class SplitExpSoftplus(nn.Module):
    def __init__(self):
        super(SplitExpSoftplus, self).__init__()
        self.softplus = nn.Softplus()

    def forward(self, x):
        mean, dispersion = torch.split(x, x.shape[1] // 2, dim=1)
        mean = CustomClamp()(mean)
        mean = Exp()(mean)
        dispersion = self.softplus(dispersion)
        return mean, dispersion


class CustomClamp(nn.Module):
    """Clamp input to [-inf, 30] and use ReLU for the gradients."""
    def __init__(self):
        super(CustomClamp, self).__init__()

    def forward(self, x):
        return 30 - torch.nn.ReLU()(30 - x)


class Exp(nn.Module):
    def __init__(self):
        super(Exp, self).__init__()

    def forward(self, x):
        return torch.exp(x)


def negative_binom_loss(gex_pred, gex_true):
    mean, dispersion = gex_pred
    eps = 1e-10  # for numerical stability

    # Clip dispersion values
    dispersion = torch.clamp(dispersion, max=1e6)

    t1 = (
            torch.lgamma(dispersion + eps)
            + torch.lgamma(gex_true + 1.0)
            - torch.lgamma(gex_true + dispersion + eps)
    )
    t2 = (dispersion + gex_true) * torch.log1p(mean / (dispersion + eps)) + (
            gex_true * (torch.log(dispersion + eps) - torch.log(mean + eps))
    )

    return torch.mean(t1 + t2)


class VanillaNN(pl.LightningModule):
    """
    A vanilla neural network for predicting ADT from GEX.
    Based on Babel, but with just GEX encoder and ADT decoder.
    """
    def __init__(self, gex_dim, adt_dim):
        super(VanillaNN, self).__init__()
        self.automatic_optimization = False

        # GEX encoder
        self.gex_encoder = nn.Sequential(
            nn.Dropout(0.6),
            nn.Linear(gex_dim, 64),
            nn.PReLU(),
            nn.Linear(64, 64),
            nn.PReLU(),
        )

        # ADT decoder
        self.adt_decoder = nn.Sequential(
            nn.Linear(64, 64),
            nn.Tanh(),
            nn.Linear(64, adt_dim),
            nn.ReLU(),
        )

    def forward(self, gex, adt):
        # encode
        gex_latent = self.gex_encoder(gex)
        # decode
        adt_from_gex = self.adt_decoder(gex_latent)
        return adt_from_gex

    def predict(self, gex):
        gex_latent = self.gex_encoder(gex)
        adt_from_gex = self.adt_decoder(gex_latent)
        return adt_from_gex

    def configure_optimizers(self):
        opt = Adam(self.parameters(), lr=0.005, weight_decay=1e-4)
        sch = torch.optim.lr_scheduler.SequentialLR(
            optimizer=opt,
            schedulers=[
                torch.optim.lr_scheduler.ConstantLR(opt, factor=1.0, total_iters=10, verbose=True),
                torch.optim.lr_scheduler.LinearLR(opt, start_factor=1.0, end_factor=1e-3,
                                                  total_iters=30, verbose=True, )
            ],
            milestones=[10],
        )
        return {"optimizer": opt, "lr_scheduler": sch, "monitor": "val_gex2adt"}

    def training_step(self, train_batch, batch_idx):
        opt = self.optimizers()
        sch = self.lr_schedulers()
        opt.optimizer.zero_grad()
        gex, adt = train_batch
        adt_from_gex = self(gex, adt)
        loss = nn.functional.mse_loss(adt_from_gex, adt)
        self.log("train_loss", loss)
        self.log("lr", opt.optimizer.param_groups[0]['lr'])
        self.manual_backward(loss)
        self.clip_gradients(opt.optimizer, 5, "norm")
        opt.step()
        if self.trainer.is_last_batch:
            sch.step()
        return loss

    def validation_step(self, val_batch, batch_idx):
        gex, adt = val_batch
        adt_from_gex = self(gex, adt)
        loss = nn.functional.mse_loss(adt_from_gex, adt)
        self.log("val_gex2adt", loss)
        return loss


class BabelDance(VanillaNN):
    """
    Babel code adapted from DANCE package: https://github.com/OmicsML/dance/blob/5edba7de34c85326bf7874cd262989f7baa2db03/examples/multi_modality/predict_modality/babel.py
    """
    def __init__(self, gex_dim, adt_dim):
        super(VanillaNN, self).__init__()

        # GEX encoder
        self.gex_encoder = nn.Sequential(
            nn.Linear(gex_dim, 64),
            nn.BatchNorm1d(64),
            nn.PReLU(),
            nn.Linear(64, 16),
            nn.BatchNorm1d(16),
            nn.PReLU(),
        )

        # GEX decoder
        self.gex_decoder = nn.Sequential(
            nn.Linear(16, 64),
            nn.BatchNorm1d(64),
            nn.PReLU(),
            nn.Linear(64, gex_dim * 2),
            # then split into mean and dispersion
            # and apply exp and softplus
            SplitExpSoftplus(),
        )

        # ADT decoder
        self.adt_decoder = nn.Sequential(
            nn.Linear(16, 64),
            nn.BatchNorm1d(64),
            nn.Tanh(),
            nn.Linear(64, adt_dim),
        )

    def forward(self, gex, adt):
        # encode
        gex_latent = self.gex_encoder(gex)
        # decode
        gex_from_gex = self.gex_decoder(gex_latent)
        adt_from_gex = self.adt_decoder(gex_latent)
        return gex_from_gex, adt_from_gex

    def configure_optimizers(self):
        opt = Adam(self.parameters(), lr=0.01)
        sch = torch.optim.lr_scheduler.ReduceLROnPlateau(opt, patience=10, factor=0.1, min_lr=1e-6, verbose=True)
        return {"optimizer": opt, "lr_scheduler": sch, "monitor": "val_gex2adt"}

    def training_step(self, train_batch, batch_idx):
        opt = self.optimizers()
        opt.optimizer.zero_grad()
        gex, adt = train_batch
        gex_from_gex, adt_from_gex = self(gex, adt)
        loss_gex2gex = negative_binom_loss(gex_from_gex, gex)
        loss_gex2adt = nn.functional.mse_loss(adt_from_gex, adt)

        self.log("gex2gex", loss_gex2gex)
        self.log("gex2adt", loss_gex2adt)
        loss = loss_gex2adt + loss_gex2gex
        self.log("train_loss", loss)
        self.log("lr", opt.optimizer.param_groups[0]['lr'])
        return loss

    def validation_step(self, val_batch, batch_idx):
        gex, adt = val_batch
        gex_from_gex, adt_from_gex = self(gex, adt)
        loss_gex2gex = negative_binom_loss(gex_from_gex, gex)
        loss_gex2adt = nn.functional.mse_loss(adt_from_gex, adt)
        loss = loss_gex2adt + loss_gex2gex
        self.log("val_loss", loss)
        self.log("val_gex2adt", loss_gex2adt)
        print("sqrt val loss", torch.sqrt(loss_gex2adt))
        return loss

